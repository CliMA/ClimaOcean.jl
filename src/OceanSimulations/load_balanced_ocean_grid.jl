using ClimaOcean.Bathymetry
using Oceananigans
using Oceananigans.Architectures: arch_array, device
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: Sizes, child_architecture
using Oceananigans.ImmersedBoundaries: immersed_cell
using KernelAbstractions: @index, @kernel
using JLD2

"""
    load_balanced_regional_grid(arch; size, 
                                latitude, 
                                longitude, 
                                z, 
                                halo, 
                                maximum_size, 
                                height_above_water, 
                                minimum_depth, 
                                interpolation_passes)

Construct a LatitudeLongitudeGrid with an ocean bathymetry interpolated from ETOPO1.
If the architecture is `Distributed` and the partition is only in one direction, the partition will be
calculated to maintain an equal number of active cells across different workers.

# Positional Arguments
======================
- `arch`: The architecture of the ocean grid.

# Keyword Arguments
===================
- `size`: The size of the grid.
- `latitude`: The latitude of the grid.
- `longitude`: The longitude of the grid.
- `z`: The z-faces of the grid.
- `halo`: The halo size of the grid. Default is `(3, 3, 3)`.
- `maximum_size`: The maximum size in the partitioned direction. In `nothing` the load balanced direction 
                  is not reduced. Default is `nothing`.
- `height_above_water`: The height above water level. Default is `1`.
- `minimum_depth`: The minimum depth of the bathymetry. Default is `10`.
- `interpolation_passes`: The number of interpolation passes. Default is `1`.
- `connected_regions_allowed`: The number of connected regions allowed in the bathymetry

# Returns
- `grid`: The load-balanced ocean grid.
"""
function load_balanced_regional_grid(arch; 
                                     size, 
                                     latitude, 
                                     longitude, 
                                     z,
                                     halo = (3, 3, 3),
                                     maximum_size = nothing,
                                     height_above_water = 1,
                                     minimum_depth = 10,
                                     connected_regions_allowed = 3, 
                                     interpolation_passes = 1,
                                     bathymetry_file = nothing)
    
    grid = LatitudeLongitudeGrid(arch;
                                 size,
                                 longitude,
                                 latitude,
                                 z,
                                 halo)

    if !isnothing(bathymetry_file)
        if isfile(bathymetry_file)
            bottom_height = jldopen(bathymetry_file)["bathymetry"]
            return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 
        else
            bottom_height = regrid_bathymetry(grid;
                                              height_above_water,
                                              minimum_depth,
                                              interpolation_passes,
                                              connected_regions_allowed)
            
            jldsave(bathymetry_file, bathymetry = Array(interior(bottom_height)))
            return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 
        end
    else
        bottom_height = regrid_bathymetry(grid;
                                          height_above_water,
                                          minimum_depth,
                                          interpolation_passes,
                                          connected_regions_allowed)
    
        return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 
    end
end

const SlabPartition = Union{Partition{<:Any, <:Nothing, <:Nothing},
                            Partition{<:Nothing, <:Any, <:Nothing}}

const SlabDistributed = Distributed{<:Any, <:Any, <:SlabPartition}

# Load balancing works only for 1D partitions!
function load_balanced_regional_grid(arch::SlabDistributed; 
                                     size, 
                                     latitude, 
                                     longitude, 
                                     z,
                                     halo = (3, 3, 3),
                                     maximum_size = nothing,
                                     height_above_water = 1,
                                     minimum_depth = 10,
                                     interpolation_passes = 1,
                                     connected_regions_allowed = 3, 
                                     bathymetry_file = nothing)
        
    child_arch = child_architecture(arch)

    # Global grid
    grid = load_balanced_regional_grid(child_arch;
                                       size,
                                       longitude,
                                       latitude,
                                       z,
                                       halo,
                                       maximum_size,
                                       height_above_water,
                                       minimum_depth,
                                       interpolation_passes,
                                       bathymetry_file,
                                       connected_regions_allowed)

    bottom_height = grid.immersed_boundary.bottom_height

    # index of the partitioned direction
    idx = arch.ranks[1] == 1 ? 2 : 1

    # Starting with the load balancing
    load_per_slab = arch_array(child_arch, zeros(Int, size[idx]))
    loop! = assess_load(device(child_arch), 512, size[idx])
    loop!(load_per_slab, grid, idx)

    load_per_slab = arch_array(CPU(), load_per_slab)
    local_N       = calculate_local_size(load_per_slab, size[idx], arch.ranks[idx])

    # Limit the maximum size such that we do not have memory issues
    redistribute_size_to_fulfill_memory_limitation!(local_N, maximum_size)

    partition = idx == 1 ? Partition(x = Sizes(local_N...)) : Partition(y = Sizes(local_N...))

    arch = Distributed(child_arch; partition)
    zonal_rank = arch.local_index[idx]

    @info "slab decomposition with " zonal_rank local_N[zonal_rank]

    grid = LatitudeLongitudeGrid(arch;
                                 size = N,
                                 longitude,
                                 latitude,
                                 z,
                                 halo)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height), active_cells_map = true)
end

@kernel function assess_load(load_per_slab, grid, idx)
    i1 = @index(Global, Linear)

    for i2 in 1:size(grid, idx)
        for k in 1:size(grid, 3)
            i = ifelse(idx == 1, (i1, i2), (i2, i1))
            @inbounds load_per_slab[i1] += ifelse(immersed_cell(i..., k, grid), 0, 1)
        end
    end
end

function calculate_local_size(load_per_slab, N, ranks)
    active_cells  = sum(load_per_slab)
    active_load   = active_cells / ranks
    local_N = zeros(Int, ranks) # fill the local N with the active load
    idx = 1
    for r in 1:ranks-1
        local_load = 0
        while local_load <= active_load
            local_load += load_per_slab[idx]
            local_N[r] += 1
            idx += 1
        end
    end

    local_N[end] = N - sum(local_N[1:end-1])

    return local_N
end

redistribute_size_to_fulfill_memory_limitation!(l, ::Nothing) = nothing

function redistribute_size_to_fulfill_memory_limitation!(l, m)
    n = length(l)
    while any(l .> m)
        x⁺, i⁺ = findmax(l)
        x⁻, i⁻ = findmin(l)
        if x⁺ == m + 1
            l[i⁺] -= 1
            l[i⁻] += 1
        else
            Δ = l[i⁺] - m
            q = Δ ÷ n
            l .+= q
            l[i⁺] -= Δ
            l[i⁻] += mod(Δ, n)
        end
    end
    while any(l .< 20)
        x⁺, i⁺ = findmax(l)
        x⁻, i⁻ = findmin(l)
        if x⁻ == 19
            l[i⁻] += 1
            l[i⁺] -= 1
        else
            Δ = 20 - l[i⁻]
            q = Δ ÷ n
            l .-= q
            l[i⁺] -= mod(Δ, n)
            l[i⁻] += Δ
        end
    end

    return nothing
end
