module Bathymetry

export regrid_bathymetry

using KernelAbstractions: @kernel, @index

using ImageMorphology

using Oceananigans
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: DistributedGrid, reconstruct_global_grid, barrier!, all_reduce
using Oceananigans.Grids: halo_size, λnodes, φnodes, x_domain, y_domain, topology
using Oceananigans.Utils: pretty_filesize, launch!
using Oceananigans.Fields: interpolate!
using Oceananigans.BoundaryConditions
using JLD2

using ClimaOcean
using OffsetArrays
using NCDatasets
using Downloads
using Printf
using Scratch

using ..DataWrangling: Metadatum, native_grid, metadata_path, download_dataset

# methods specific to bathymetric datasets are added within dataset modules

"""
    regrid_bathymetry(target_grid;
                      height_above_water = nothing,
                      minimum_depth = 0,
                      dir = download_bathymetry_cache,
                      url = ClimaOcean.Bathymetry.etopo_url,
                      filename = "ETOPO_2022_v1_60s_N90W180_surface.nc",
                      interpolation_passes = 1,
                      major_basins = 1)

Return bathymetry associated with the NetCDF file at `path = joinpath(dir, filename)` regridded onto `target_grid`.
If `path` does not exist, then a download is attempted from `joinpath(url, filename)`.

Arguments
=========

- `target_grid`: grid to interpolate the bathymetry onto.

Keyword Arguments
=================

- `height_above_water`: limits the maximum height of above-water topography (where ``h > 0``) before interpolating.
                        Default: `nothing`, which implies that the original topography is retained.

- `minimum_depth`: minimum depth for the shallow regions, defined as a positive value.
                   `h > - minimum_depth` is considered land. Default: 0.

- `dir`: directory of the bathymetry-containing file. Default: `download_bathymetry_cache`.

- `filename`: file containing bathymetric data. Must be netCDF with fields:
  1. `lat` vector of latitude nodes
  2. `lon` vector of longitude nodes
  3. `z` matrix of depth values

- `interpolation_passes`: regridding/interpolation passes. The bathymetry is interpolated in
                          `interpolation_passes - 1` intermediate steps. The more the interpolation
                          steps the smoother the final bathymetry becomes.

  Example
  =======

  Interpolating from a 400 x 200 grid to a 100 x 100 grid in 4 passes involves:

  * 400 x 200 → 325 x 175
  * 325 x 175 → 250 x 150
  * 250 x 150 → 175 x 125
  * 175 x 125 → 100 x 100

  If _coarsening_ the original grid, linear interpolation in passes is equivalent to
  applying a smoothing filter, with more passes increasing the strength of the filter.
  If _refining_ the original grid, additional passes do not help and no intermediate
  steps are performed.

- `major_basins`: Number of "independent major basins", or fluid regions fully encompassed by land,
                  that are retained by [`remove_minor_basins!`](@ref). Basins are removed by order of size:
                  the smallest basins are removed first. `major_basins = 1` retains only the largest basin.
                  If `Inf` then no basins are removed. Default: 1.
"""
function regrid_bathymetry(target_grid, metadata;
                           height_above_water = nothing,
                           minimum_depth = 0,
                           interpolation_passes = 1,
                           major_basins = 1) # Allow an `Inf` number of "lakes"

    download_dataset(metadata)

    if isinteger(interpolation_passes)
        interpolation_passes = convert(Int, interpolation_passes)
    end

    if interpolation_passes isa Nothing || !isa(interpolation_passes, Int) || interpolation_passes ≤ 0
        return throw(ArgumentError("interpolation_passes has to be an integer ≥ 1"))
    end

    arch = architecture(target_grid)

    bathymetry_native_grid = native_grid(metadata, arch; halo = (10, 10, 1))
    FT = eltype(target_grid)

    filepath = metadata_path(metadata)
    dataset = Dataset(filepath, "r")

    z_data = convert(Array{FT}, dataset["z"][:, :])
    close(dataset)

    if !isnothing(height_above_water)
        # Overwrite the height of cells above water.
        # This has an impact on reconstruction. Greater height_above_water reduces total
        # wet area by biasing coastal regions to land during bathymetry regridding.
        land = z_data .> 0
        z_data[land] .= height_above_water
    end

    native_z = Field{Center, Center, Nothing}(bathymetry_native_grid)
    set!(native_z, z_data)
    fill_halo_regions!(native_z)


    target_z = interpolate_bathymetry_in_passes(native_z, target_grid;
                                                passes = interpolation_passes)

    if minimum_depth > 0
        launch!(arch, target_grid, :xy, _enforce_minimum_depth!, target_z, minimum_depth)
    end

    if major_basins < Inf
        remove_minor_basins!(target_z, major_basins)
    end

    fill_halo_regions!(target_z)

    return target_z
end

function regrid_bathymetry(target_grid; dataset = ClimaOcean.ETOPO.ETOPO2022(), kw...)
    metadatum = Metadatum(:bottom_height; dataset)

    return regrid_bathymetry(target_grid, metadatum; kw...)
end

@kernel function _enforce_minimum_depth!(target_z, minimum_depth)
    i, j = @index(Global, NTuple)
    z = @inbounds target_z[i, j, 1]

    # Fix active cells to be at least `-minimum_depth`.
    active = z < 0 # it's a wet cell
    z = ifelse(active, min(z, -minimum_depth), z)

    @inbounds target_z[i, j, 1] = z
end

# Here we can either use `regrid!` (three dimensional version) or `interpolate!`.
function interpolate_bathymetry_in_passes(native_z, target_grid;
                                          passes = 10)

    gridtype = target_grid isa TripolarGrid ? "TripolarGrid" :
               target_grid isa LatitudeLongitudeGrid ? "LatitudeLongitudeGrid" :
               target_grid isa RectilinearGrid ? "RectilinearGrid" :
               error("unknown target grid type")

    Nλt, Nφt = Nt = size(target_grid)
    Nλn, Nφn = Nn = size(native_z)

    resxt = minimum_xspacing(target_grid)
    resyt = minimum_yspacing(target_grid)

    resxn = minimum_xspacing(native_z.grid)
    resyn = minimum_yspacing(native_z.grid)

    # Check whether we are refining the grid in any directions.
    # If so, skip interpolation passes, as they are not needed.
    if resxt < resxn || resyt < resyn
        target_z = Field{Center, Center, Nothing}(target_grid)

        @info "Interpolating bathymetry of size $Nn onto a $gridtype target grid of size $Nt"
        interpolate!(target_z, native_z)

        if passes > 0
            passes_str = passes == 1 ? "1 interpolation pass" : "$passes interpolation passes"

            @info string("Skipping extra $passes_str for bathymetry.", '\n',
                        "Extra interpolation passes may only be used to coarsen bathymetry", '\n',
                        "and require that the bathymetry is finer than the target grid in", '\n',
                        "both horizontal directions.", '\n',
                        "Grid info:", '\n',
                        "target grid", '\n',
                        "├── minimum x-spacings ", @sprintf("%.3e", resxt), " m", '\n',
                        "└── minimum y-spacings ", @sprintf("%.3e", resyt), " m", '\n',
                        "bathymetry grid", '\n',
                        "├── minimum x-spacings ", @sprintf("%.3e", resxn), " m", '\n',
                        "└── minimum y-spacings ", @sprintf("%.3e", resyn), " m")
        end
        return target_z
    end

    # Interpolate in passes
    latitude  = y_domain(native_z.grid)
    longitude = x_domain(native_z.grid)

    ΔNλ = floor((Nλn - Nλt) / passes)
    ΔNφ = floor((Nφn - Nφt) / passes)

    Nλ = [Nλn - ΔNλ * pass for pass in 1:passes-1]
    Nφ = [Nφn - ΔNφ * pass for pass in 1:passes-1]

    Nλ = Int[Nλ..., Nλt]
    Nφ = Int[Nφ..., Nφt]

    old_z  = native_z
    TX, TY = topology(target_grid)

    @info "Interpolation passes of bathymetry size $(size(old_z)) onto a $gridtype target grid of size $Nt:"
    for pass = 1:passes - 1
        new_size = (Nλ[pass], Nφ[pass], 1)
        @info "    pass $pass to size $new_size"

        new_grid = LatitudeLongitudeGrid(architecture(target_grid), Float32,
                                         size = new_size,
                                         latitude = (latitude[1],  latitude[2]),
                                         longitude = (longitude[1], longitude[2]),
                                         z = (0, 1),
                                         topology = (TX, TY, Bounded))

        new_z = Field{Center, Center, Nothing}(new_grid)

        interpolate!(new_z, old_z)
        old_z = new_z
    end

    new_size = (Nλ[passes], Nφ[passes], 1)
    @info "    pass $passes to size $new_size"
    target_z = Field{Center, Center, Nothing}(target_grid)
    interpolate!(target_z, old_z)

    return target_z
end

# Regridding bathymetry for distributed grids, we handle the whole process
# on just one rank, and share the results with the other processors.
function regrid_bathymetry(target_grid::DistributedGrid, args...; kw...)
    global_grid = reconstruct_global_grid(target_grid)
    global_grid = on_architecture(CPU(), global_grid)
    arch = architecture(target_grid)
    Nx, Ny, _ = size(global_grid)

    # If all ranks open a gigantic bathymetry and the memory is
    # shared, we could easily have OOM errors.
    # We perform the reconstruction only on rank 0 and share the result.
    bottom_height = if arch.local_rank == 0
        bottom_field = regrid_bathymetry(global_grid, args...; kw...)
        bottom_field.data[1:Nx, 1:Ny, 1]
    else
        zeros(Nx, Ny)
    end

    # Synchronize
    Oceananigans.DistributedComputations.global_barrier(arch.communicator)

    # Share the result (can we share SubArrays?)
    bottom_height = all_reduce(+, bottom_height, arch)

    # Partition the result
    local_bottom_height = Field{Center, Center, Nothing}(target_grid)
    set!(local_bottom_height, bottom_height)
    fill_halo_regions!(local_bottom_height)

    return local_bottom_height
end


"""
    remove_minor_basins!(z_data, keep_major_basins)

Remove independent basins from the bathymetry data stored in `z_data` by identifying connected regions
below sea level. Basins are removed from smallest to largest until only `keep_major_basins` remain.

Arguments
=========

- `z_data`: A 2D array representing the bathymetry data.
- `keep_major_basins`: The maximum number of connected regions to keep.
                       If `Inf` is provided then all connected regions are kept.

"""
function remove_minor_basins!(zb::Field, keep_major_basins)
    zb_cpu = on_architecture(CPU(), zb)
    TX     = topology(zb_cpu.grid, 1)

    Nx, Ny, _ = size(zb_cpu.grid)
    zb_data   = maybe_extend_longitude(zb_cpu, TX()) # Outputs a 2D AbstractArray

    remove_minor_basins!(zb_data, keep_major_basins)
    set!(zb, zb_data[1:Nx, 1:Ny])

    return zb
end

maybe_extend_longitude(zb_cpu, tx) = interior(zb_cpu, :, :, 1)

# Since the strel algorithm in `remove_major_basins` does not recognize periodic boundaries,
# before removing connected regions, we extend the longitude direction if it is periodic.
# An extension of half the domain is enough.
function maybe_extend_longitude(zb_cpu, ::Periodic)
    Nx = size(zb_cpu, 1)
    nx = Nx ÷ 2

    zb_data   = zb_cpu.data[1:Nx, :, 1]
    zb_parent = zb_data.parent

    # Add information on the LHS and to the RHS
    zb_parent = vcat(zb_parent[nx:Nx, :], zb_parent, zb_parent[1:nx, :])

    # Update offsets
    yoffsets = zb_cpu.data.offsets[2]
    xoffsets = - nx

    return OffsetArray(zb_parent, xoffsets, yoffsets)
end

remove_major_basins!(zb::OffsetArray, keep_minor_basins) =
    remove_minor_basins!(zb.parent, keep_minor_basins)

function remove_minor_basins!(zb, keep_major_basins)

    if !isfinite(keep_major_basins)
        throw(ArgumentError("`keep_major_basins` must be a finite number!"))
    end

    if keep_major_basins < 1
        throw(ArgumentError("keep_major_basins must be larger than 0."))
    end

    water = zb .< 0

    connectivity = ImageMorphology.strel(water)
    labels = ImageMorphology.label_components(connectivity)

    total_elements = zeros(maximum(labels))
    label_elements = zeros(maximum(labels))

    for e in 1:lastindex(total_elements)
        total_elements[e] = length(labels[labels .== e])
        label_elements[e] = e
    end

    mm_basins = [] # major basins indexes
    m = 1

    # We add basin indexes until we reach the specified number (m == keep_major_basins) or
    # we run out of basins to keep -> isempty(total_elements)
    while (m <= keep_major_basins) && !isempty(total_elements)
        next_maximum = findfirst(x -> x == maximum(total_elements), total_elements)
        push!(mm_basins, label_elements[next_maximum])
        deleteat!(total_elements, next_maximum)
        deleteat!(label_elements, next_maximum)
        m += 1
    end

    labels = map(Float64, labels)

    for ℓ = 1:maximum(labels)
        remove_basin = all(ℓ != m for m in mm_basins)
        if remove_basin
            labels[labels .== ℓ] .= NaN # Regions to remove
        end
    end

    # Flatten minor basins, corresponding to regions where `labels == NaN`
    zb[isnan.(labels)] .= 0

    return nothing
end

end # module
