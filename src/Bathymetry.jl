module Bathymetry

export regrid_bathymetry, retrieve_bathymetry

using ImageMorphology
using ..DataWrangling: download_progress

using Oceananigans
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.Grids: halo_size, λnodes, φnodes
using Oceananigans.Grids: x_domain, y_domain
using Oceananigans.Grids: topology
using Oceananigans.Utils: pretty_filesize, launch!
using Oceananigans.Fields: interpolate!
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index
using JLD2

using NCDatasets
using Downloads
using Printf
using Scratch

download_bathymetry_cache::String = ""
function __init__()
    global download_bathymetry_cache = @get_scratch!("Bathymetry")
end

"""
    regrid_bathymetry(target_grid;
                      url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf", 
                      height_above_water = <none>,
                      minimum_depth = 0,
                      dir = download_cache,
                      filename = "ETOPO_2022_v1_60s_N90W180_surface.nc")

Regrid bathymetry associated with the NetCDF file at `path = joinpath(dir, filename)` to `target_grid`.
If `path` does not exist, then a download is attempted from `joinpath(url, filename)`.

Arguments:
==========

- target_grid: grid to interpolate onto

Keyword Arguments:
==================

- height_above_water: limits the maximum height of above-water topography (where h > 0). If
                      `nothing` the original topography is retained

- minimum_depth: minimum depth for the shallow regions, defined as a positive value. 
                 `h > - minimum_depth` will be considered land

- dir: directory of the bathymetry-containing file

- filename: file containing bathymetric data. Must be netcdf with fields:
            (1) `lat` vector of latitude nodes
            (2) `lon` vector of longitude nodes
            (3) `z` matrix of depth values

- interpolation_passes: regridding/interpolation passes. The bathymetry is interpolated in
                        `interpolation_passes - 1` intermediate steps. With more steps the 
                        final bathymetry will be smoother.
                        Example: interpolating from a 400x200 grid to a 100x100 grid in 4 passes will involve
                        - 400x200 -> 325x175
                        - 325x175 -> 250x150
                        - 250x150 -> 175x125
                        - 175x125 -> 100x100
                        If _coarsening_ the original grid, linear interpolation in passes is equivalent to 
                        applying a smoothing filter, with more passes increasing the strength of the filter.
                        If _refining_ the original grid, additional passes will not help and no intermediate
                        steps will be performed.

- major_basins: Number of "independent major basins", or fluid regions fully encompassed by land,
                that are retained by [`remove_minor_basins!`](@ref). Basins are removed by order of size:
                the smallest basins are removed first. `major_basins=1` will retain only the largest basin.
                Default: `Inf`, which does not remove any basins.
"""
function regrid_bathymetry(target_grid;
                           height_above_water = nothing,
                           minimum_depth = 0,
                           dir = download_bathymetry_cache,
                           url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf", 
                           filename = "ETOPO_2022_v1_60s_N90W180_surface.nc",
                           interpolation_passes = 1,
                           major_basins = Inf) # Allow an `Inf` number of ``lakes''

    filepath = joinpath(dir, filename)

    if isfile(filepath)
        @info "Regridding bathymetry from existing file $filepath."
    else
        @info "Downloading bathymetry..."
        if !ispath(dir)
            @info "Making bathymetry directory $dir..."
            mkdir(dir)
        end

        fileurl = joinpath(url, filename)

        try 
            Downloads.download(fileurl, filepath; progress=download_progress, verbose=true)
        catch 
            cmd = `wget --no-check-certificate -O $filepath $fileurl`
            run(cmd)
        end
    end

    dataset = Dataset(filepath)

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
    arch = child_architecture(architecture(target_grid))
    φ₁, φ₂ = y_domain(target_grid)
    λ₁, λ₂ = x_domain(target_grid)

    if λ₁ < 0 || λ₂ > 360
        throw(ArgumentError("Cannot regrid bathymetry between λ₁ = $(λ₁) and λ₂ = $(λ₂).
                             Bathymetry data is defined on longitudes spanning λ = (0, 360)."))
    end

    # Calculate limiting indices on the bathymetry grid
    i₁ = searchsortedfirst(λ_data, λ₁)
    i₂ = searchsortedfirst(λ_data, λ₂) - 1
    ii = i₁:i₂

    j₁ = searchsortedfirst(φ_data, φ₁)
    j₂ = searchsortedfirst(φ_data, φ₂) - 1
    jj = j₁:j₂

    # Restrict bathymetry _data to region of interest
    λ_data = λ_data[ii]
    φ_data = φ_data[jj]
    z_data = z_data[ii, jj]

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

    λ₁_data = λ_data[1]   - Δλ / 2
    λ₂_data = λ_data[end] + Δλ / 2
    φ₁_data = φ_data[1]   - Δφ / 2
    φ₂_data = φ_data[end] + Δφ / 2

    Nxn = length(λ_data)
    Nyn = length(φ_data)
    Nzn = 1

    native_grid = LatitudeLongitudeGrid(arch;
                                        size = (Nxn, Nyn, Nzn),
                                        latitude  = (φ₁_data, φ₂_data),
                                        longitude = (λ₁_data, λ₂_data),
                                        z = (0, 1),
                                        halo = (10, 10, 1))

    native_z = Field{Center, Center, Nothing}(native_grid)
    set!(native_z, z_data)

    target_z = interpolate_bathymetry_in_passes(native_z, target_grid; 
                                                passes = interpolation_passes,
                                                minimum_depth)

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
function interpolate_bathymetry_in_passes(native_z, target_grid; 
                                          passes = 10,
                                          minimum_depth = 0)
    Nλt, Nφt = Nt = size(target_grid)
    Nλn, Nφn = Nn = size(native_z)

    # Check whether we are coarsening the grid in any directions.
    # If so, skip interpolation passes.
    if Nλt > Nλn || Nφt > Nφn
        target_z = Field{Center, Center, Nothing}(target_grid)
        interpolate!(target_z, native_z)
        @info string("Skipping passes for interplating bathymetry of size $Nn", '\n',
                     "to target grid of size $Nt. Interpolation passes may only", '\n',
                     "be used to refine bathymetryand requires that the bathymetry", '\n',
                     "is larger than the target grid in both horizontal directions.")
        return target_z
    end
 
    # Interpolate in passes
    latitude  = y_domain(target_grid)
    longitude = x_domain(target_grid)

    ΔNλ = floor((Nλn - Nλt) / passes)
    ΔNφ = floor((Nφn - Nφt) / passes)

    Nλ = [Nλn - ΔNλ * pass for pass in 1:passes-1]
    Nφ = [Nφn - ΔNφ * pass for pass in 1:passes-1]

    Nλ = Int[Nλ..., Nλt]
    Nφ = Int[Nφ..., Nφt]

    old_z  = native_z
    TX, TY = topology(target_grid)

    for pass = 1:passes - 1
        new_size = (Nλ[pass], Nφ[pass], 1)

        @debug "Bathymetry interpolation pass $pass with size $new_size"

        new_grid = LatitudeLongitudeGrid(architecture(target_grid),
                                         size = new_size, 
                                         latitude = (latitude[1],  latitude[2]), 
                                         longitude = (longitude[1], longitude[2]), 
                                         z = (0, 1),
                                         topology = (TX, TY, Bounded))

        new_z = Field{Center, Center, Nothing}(new_grid)

        interpolate!(new_z, old_z)
        old_z = new_z
    end

    target_z = Field{Center, Center, Nothing}(target_grid)
    interpolate!(target_z, old_z)

    return target_z
end

"""
    remove_minor_basins!(z_data, keep_major_basins)

Remove independent basins from the bathymetry data stored in `z_data` by identifying connected regions
below sea level. Basins are removed from smallest to largest until only `keep_major_basins` remain.

Arguments
=========

- `z_data`: A 2D array representing the bathymetry data.
- `keep_major_basins`: The maximum number of connected regions to keep. 
                       Default is `Inf`, which means all connected regions are kept.

"""
function remove_minor_basins!(Z::Field, keep_major_basins)
    Zi = interior(Z, :, :, 1)
    Zi_cpu = on_architecture(CPU(), Zi)
    remove_minor_basins!(Zi_cpu, keep_major_basins)
    set!(Z, Zi_cpu)

    return Z
end

function remove_minor_basins!(Z, keep_major_basins)

    if !isfinite(keep_major_basins)
        throw(ArgumentError("`keep_major_basins` must be a finite number!"))
    end

    if keep_major_basins < 1
        throw(ArgumentError("keep_major_basins must be larger than 0."))
    end

    water = Z .< 0
    
    connectivity = ImageMorphology.strel(water)
    labels = ImageMorphology.label_components(connectivity)
    
    total_elements = zeros(maximum(labels))
    label_elements = zeros(maximum(labels))

    for e in 1:lastindex(total_elements)
        total_elements[e] = length(labels[labels .== e])
        label_elements[e] = e
    end
        
    mm_basins = [] # major basins indexes
    for m = 1:keep_major_basins
        next_maximum = findfirst(x -> x == maximum(total_elements), total_elements)
        push!(mm_basins, label_elements[next_maximum])
        total_elements = filter(x -> x != total_elements[next_maximum], total_elements)
        label_elements = filter(x -> x != label_elements[next_maximum], label_elements)
    end
        
    labels = map(Float64, labels)

    for ℓ = 1:maximum(labels)
        remove_basin = all(ℓ != m for m in mm_basins)
        if remove_basin
            labels[labels .== ℓ] .= 1e10 # large number
        end
    end

    # Flatten land
    labels[labels .<  1e10] .= 0 
    labels[labels .== 1e10] .= NaN

    Z .+= labels
    Z[isnan.(Z)] .= 0

    return nothing
end

"""
    retrieve_bathymetry(grid, filename; kw...)

Retrieve the bathymetry data from a file or generate it using a grid and save it to a file.

# Arguments
============

- `grid`: The grid used to generate the bathymetry data.
- `filename`: The name of the file to read or save the bathymetry data.
- `kw...`: Additional keyword arguments.

# Returns
===========
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

