using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: Sizes
using Oceananigans.ImmersedBoundaries: immersed_cell
using KernelAbstractions: @index, @kernel

function LoadBalancedOceanGrid(arch; 
                               size, 
                               latitude, 
                               longitude, 
                               z,
                               halo = (3, 3, 3),
                               kw...)

    child_arch = child_architecture(arch)
    
    grid = LatitudeLongitudeGrid(child_arch;
                                 size,
                                 longitude,
                                 latitude,
                                 z,
                                 halo)

    bottom_height = regrid_bathymetry(grid, 
                                      height_above_water = 1,
                                      minimum_depth = 10,
                                      interpolation_passes = 25)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 
end

function LoadBalancedOceanGrid(arch::Distributed; 
                               size, 
                               latitude, 
                               longitude, 
                               z,
                               halo = (3, 3, 3),
                               maximum_Nx = 1150)

    child_arch = child_architecture(arch)

    # Global grid
    grid = LatitudeLongitudeGrid(child_arch;
                                 size,
                                 longitude,
                                 latitude,
                                 z,
                                 halo)

    bottom_height = regrid_bathymetry(grid, 
                                      height_above_water = 1,
                                      minimum_depth = 10,
                                      interpolation_passes = 25)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

    # Starting with the load balancing
    load_per_x_slab = arch_array(child_arch, zeros(Int, size[1]))
    loop! = assess_x_load(device(child_arch), 512, size[1])
    loop!(load_per_x_slab, grid)

    load_per_x_slab = arch_array(CPU(), load_per_x_slab)
    local_Nx        = calculate_local_size(load_per_x_slab, N[1], arch.ranks[1])

    # We cannot have Nx > 650 if Nranks = 32 otherwise we incur in memory limitations,
    # so for a small number of GPUs we are limited in the load balancing
    redistribute_size_to_fulfill_memory_limitation!(local_Nx, maximum_Nx)

    arch = Distributed(child_arch, partition = Partition(x = Sizes(local_Nx...)))
    zonal_rank = arch.local_index[1]

    @info "slab decomposition with " zonal_rank local_Nx[zonal_rank], arch

    @show underlying_grid = LatitudeLongitudeGrid(arch;
                                                  size = N,
                                                  longitude = (-180, 180),
                                                  latitude = latitude,
                                                  halo = (7, 7, 7),
                                                  z = z_faces)

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height), active_cells_map = true)
end

@kernel function assess_x_load(load_per_slab, ibg)
    i = @index(Global, Linear)

    @unroll for j in 1:size(ibg, 2)
        @unroll for k in 1:size(ibg, 3)
            @inbounds load_per_slab[i] += ifelse(immersed_cell(i, j, k, ibg), 0, 1)
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

function redistribute_size_to_fulfill_memory_limitation!(l, m = 700)
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
end
