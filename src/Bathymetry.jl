module Bathymetry

export regrid_bathymetry, retrieve_bathymetry

using ImageMorphology
using ..DataWrangling: download_progress

using Oceananigans
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.DistributedComputations: DistributedGrid, reconstruct_global_grid, barrier!, all_reduce
using Oceananigans.Grids: halo_size, λnodes, φnodes
using Oceananigans.Grids: x_domain, y_domain
using Oceananigans.Grids: topology
using Oceananigans.Utils: pretty_filesize, launch!
using Oceananigans.Fields: interpolate!
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index
using JLD2

using OffsetArrays
using ClimaOcean

using NCDatasets
using Downloads
using Printf
using Scratch

download_bathymetry_cache::String = ""
function __init__()
    global download_bathymetry_cache = @get_scratch!("Bathymetry")
end

# etopo_url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf/" *
#              "ETOPO_2022_v1_60s_N90W180_surface.nc"

etopo_url = "https://www.dropbox.com/scl/fi/6pwalcuuzgtpanysn4h6f/" *
            "ETOPO_2022_v1_60s_N90W180_surface.nc?rlkey=2t7890ruyk4nd5t5eov5768lt&st=yfxsy1lu&dl=0"

struct InterpolationPasses
    passes :: Int
end

interpolating(interp::InterpolationPasses, pass, z₁, grid₁, z₀, grid₀) = pass > interp.passes
interpolating(criteria::Tuple, args...) = any(interpolating(crit, args...) for crit in criteria)

struct MustRefine end

function interpolating(::MustRefine, p, z₁, grid₁, z₀, grid₀)
    both_x_regular = grid₁ isa XRegularLLG && grid₀ isa XRegularLLG
    both_y_regular = grid₁ isa YRegularLLG && grid₀ isa YRegularLLG

    min_Δx₁ = minimum_xspacing(grid₁)
    min_Δx₀ = minimum_xspacing(grid₀)
    min_Δy₁ = minimum_yspacing(grid₁)
    min_Δy₀ = minimum_yspacing(grid₀)

    refining_in_x = both_x_regular && min_Δx₁ <= min_Δx₀
    refining_in_y = both_y_regular && min_Δy₁ <= min_Δy₀

    return refining_in_x || refining_in_y
end

tupleit(a) = tuple(a)
tupleit(a::Tuple) = a
tupleit(a::Vector) = tuple(a...)

"""
    regrid_bathymetry(target_grid;
                      height_above_water = nothing,
                      minimum_depth = 0,
                      dir = download_bathymetry_cache,
                      url = ClimaOcean.Bathymetry.etopo_url,
                      filename = "ETOPO_2022_v1_60s_N90W180_surface.nc",
                      regridding_criteria = default_criteria(grid),
                      major_basins = 1)

Return bathymetry associated with the NetCDF file at `path = joinpath(dir, filename)` regridded onto `target_grid`.
If `path` does not exist, then a download is attempted from `joinpath(url, filename)`.

Arguments
=========

- `target_grid`: grid to interpolate the bathymetry onto.

Keyword Arguments
=================

- `height_above_water`: limits the maximum height of above-water topography (where ``h > 0``) before inetrpolating.
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
function regrid_bathymetry(target_grid::AbstractGrid;
                           height_above_water = nothing,
                           minimum_depth = 0,
                           dir = download_bathymetry_cache,
                           url = etopo_url,
                           filename = "ETOPO_2022_v1_60s_N90W180_surface.nc",
                           regridding_criteria = InterpolationPasses(1),
                           major_basins = 1) # Allow an `Inf` number of "lakes"

    filepath = joinpath(dir, filename)

    # No need for @root here, because only rank 0 accesses this function
    if !isfile(filepath)
        Downloads.download(url, filepath; progress=download_progress)
    end

    dataset = Dataset(filepath, "r")

    FT = eltype(target_grid)

    φ_data = dataset["lat"][:]
    λ_data = dataset["lon"][:]
    z_data = convert(Array{FT}, dataset["z"][:, :])

    # Convert longitude from (-180, 180) to (0, 360)
    λ_data .+= 180
    Nhx    = size(z_data, 1)
    z_data = circshift(z_data, (Nhx ÷ 2, 0))

    close(dataset)

    # Diagnose target grid information
    arch   = architecture(target_grid)
    λ_data = λ_data |> Array{BigFloat}
    φ_data = φ_data |> Array{BigFloat}
    
    if !isnothing(height_above_water)
        # Overwrite the height of cells above water.
        # This has an impact on reconstruction. Greater height_above_water reduces total
        # wet area by biasing coastal regions to land during bathymetry regridding.
        land = z_data .> 0
        z_data[land] .= height_above_water
    end

    # Infer the "native grid" of the bathymetry data and make a bathymetry field.
    Δλ = λ_data[2] - λ_data[1]
    Δφ = φ_data[2] - φ_data[1]

    λ₁_data = convert(Float64, λ_data[1]   - Δλ / 2)
    λ₂_data = convert(Float64, λ_data[end] + Δλ / 2)
    φ₁_data = convert(Float64, φ_data[1]   - Δφ / 2)
    φ₂_data = convert(Float64, φ_data[end] + Δφ / 2)

    Nxn = length(λ_data)
    Nyn = length(φ_data)
    Nzn = 1

    native_grid = LatitudeLongitudeGrid(arch, Float32;
                                        size = (Nxn, Nyn, Nzn),
                                        latitude  = (φ₁_data, φ₂_data),
                                        longitude = (λ₁_data, λ₂_data),
                                        z = (0, 1),
                                        halo = (10, 10, 1))

    native_z = Field{Center, Center, Nothing}(native_grid)
    set!(native_z, z_data)
    fill_halo_regions!(native_z)

    if regridding_criteria isa Number
        regridding_criteria = tuple(InterpolationPasses(regridding_criteria))
    else
        regridding_criteria = tupleit(regridding_criteria)
    end
    regridding_criteria = tuple(MustRefine(), regridding_criteria...)
    target_z = regrid_bathymetry(native_z, target_grid, regridding_criteria)

    if minimum_depth > 0
        zi = interior(target_z, :, :, 1)

        # Set the height of cells with z > -mininum_depth to z=0.
        # (In-place + GPU-friendly)
        zi .*= zi .<= - minimum_depth
    end

    if major_basins < Inf
        remove_minor_basins!(target_z, major_basins)
    end

    fill_halo_regions!(target_z)

    return target_z
end

# Here we can either use `regrid!` (three dimensional version) or `interpolate`
function regrid_bathymetry(native_z::Field, target_grid::AbstractGrid, regridding_criteria)
                                          
    Nλt, Nφt = Nt = size(target_grid)
    Nλn, Nφn = Nn = size(native_z)

    # Interpolate in passes
    native_grid = native_z.grid
    latitude  = y_domain(native_z.grid)
    longitude = x_domain(native_z.grid)

    ΔNλ = floor((Nλn - Nλt) / passes)
    ΔNφ = floor((Nφn - Nφt) / passes)

    Nλ = [Nλn - ΔNλ * pass for pass in 1:passes-1]
    Nφ = [Nφn - ΔNφ * pass for pass in 1:passes-1]

    Nλ = Int[Nλ..., Nλt]
    Nφ = Int[Nφ..., Nφt]

    TX, TY = topology(target_grid)
    target_z = native_z
    pass = 1

    while interpolating(regridding_criteria, pass, target_z, target_grid, native_z, native_grid)
        new_size = (Nλ[pass], Nφ[pass], 1)
        @debug "Bathymetry interpolation pass $pass with size $new_size"
        new_grid = LatitudeLongitudeGrid(architecture(target_grid), Float32,
                                         size = new_size,
                                         latitude = (latitude[1],  latitude[2]),
                                         longitude = (longitude[1], longitude[2]),
                                         z = (0, 1),
                                         topology = (TX, TY, Bounded))

        old_z = target_z
        target_z = Field{Center, Center, Nothing}(new_grid)
        interpolate!(target_z, old_z)
        pass += 1
    end

    return target_z
end

# Regridding bathymetry for distributed grids, we handle the whole process
# on just one rank, and share the results with the other processors.
function regrid_bathymetry(target_grid::DistributedGrid; kw...)
    global_grid = reconstruct_global_grid(target_grid)
    global_grid = on_architecture(CPU(), global_grid)
    arch = architecture(target_grid)
    Nx, Ny, _ = size(global_grid)

    # If all ranks open a gigantic bathymetry and the memory is 
    # shared, we could easily have OOM errors. 
    # We perform the reconstruction only on rank 0 and share the result.
    bottom_height = if arch.local_rank == 0
        bottom_field = regrid_bathymetry(global_grid; kw...)
        bottom_field.data[1:Nx, 1:Ny, 1]
    else
        zeros(Nx, Ny)
    end

    # Synchronize
    ClimaOcean.global_barrier(arch.communicator)

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

"""
    retrieve_bathymetry(grid, filename; kw...)

Retrieve the bathymetry data from a file or generate it using a grid and save it to a file.

Arguments
=========

- `grid`: The grid used to generate the bathymetry data.
- `filename`: The name of the file to read or save the bathymetry data.
- `kw...`: Additional keyword arguments.

Returns
=======

- `bottom_height`: The retrieved or generated bathymetry data.

If the specified file exists, the function reads the bathymetry data from the file. 
Otherwise, it generates the bathymetry data using the provided grid and saves it to the file before returning it.
"""
function retrieve_bathymetry(grid, filename; kw...) 
    
    if isfile(filename)
        bottom_height = jldopen(filename)["bathymetry"]
    else
        bottom_height = regrid_bathymetry(grid; kw...)
        jldsave(filename, bathymetry = Array(interior(bottom_height)))
    end

    return bottom_height
end

retrieve_bathymetry(grid, ::Nothing; kw...) = regrid_bathymetry(grid; kw...)
retrieve_bathymetry(grid; kw...)            = regrid_bathymetry(grid; kw...)

end # module
