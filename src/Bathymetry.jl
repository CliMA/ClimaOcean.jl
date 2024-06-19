module Bathymetry

export regrid_bathymetry, retrieve_bathymetry

using ImageMorphology
using ..DataWrangling: download_progress

using Oceananigans
using Oceananigans.Architectures: architecture
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

"""
    regrid_bathymetry(target_grid;
                      url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf", 
                      height_above_water = <none>,
                      minimum_depth = 0,
                      dir = joinpath(@__DIR__, "..", "data"),
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

- connected_regions_allowed: number of ``connected regions'' allowed in the bathymetry. Connected regions are fluid 
                             regions that are fully encompassed by land (for example the ocean is one connected region).
                             Default is `Inf`. If a value < `Inf` is specified, connected regions will be preserved in order
                             of how many active cells they contain.
"""
function regrid_bathymetry(target_grid;
                           height_above_water = nothing,
                           minimum_depth = 0,
                           dir = joinpath(@__DIR__, "..", "data"),
                           url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf", 
                           filename = "ETOPO_2022_v1_60s_N90W180_surface.nc",
                           interpolation_passes = 1,
                           connected_regions_allowed = Inf) # Allow an `Inf` number of ``lakes''

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
        Downloads.download(fileurl, filepath; progress=download_progress, verbose=true)
    end

    dataset = Dataset(filepath)

    FT = eltype(target_grid)

    φ_data = dataset["lat"][:]
    λ_data = dataset["lon"][:]
    h_data = convert.(FT, dataset["z"][:, :])

    # Convert longitude to 0 - 360?
    λ_data .+= 180
    nhx      = size(h_data, 1)
    h_data   = circshift(h_data, (nhx ÷ 2, 1))

    close(dataset)

    # Diagnose target grid information
    arch = architecture(target_grid)
    φ₁, φ₂ = y_domain(target_grid)
    λ₁, λ₂ = x_domain(target_grid)

    # Calculate limiting indices on the bathymetry grid
    i₁ = searchsortedfirst(λ_data, λ₁)
    i₂ = searchsortedfirst(λ_data, λ₂) - 1
    ii = i₁:i₂

    j₁ = searchsortedfirst(φ_data, φ₁)
    j₂ = searchsortedfirst(φ_data, φ₂) - 1
    jj = j₁:j₂

    # Remark if `target_grid` is not perfectly nested within the bathymetry grid
    Δλ = λ_data[2] - λ_data[1]
    Δφ = φ_data[2] - φ_data[1]

    λ₁_data = λ_data[i₁] - Δλ / 2
    λ₂_data = λ_data[i₂] + Δλ / 2
    φ₁_data = φ_data[j₁] - Δφ / 2
    φ₂_data = φ_data[j₂] + Δφ / 2

    λ₁ ≈ λ₁_data || @warn "The westernmost meridian of `target_grid` $λ₁ does not coincide with " *
                          "the closest meridian of the bathymetry grid, $λ₁_data."
    λ₂ ≈ λ₂_data || @warn "The easternmost meridian of `target_grid` $λ₂ does not coincide with " *
                          "the closest meridian of the bathymetry grid, $λ₂_data."
    φ₁ ≈ φ₁_data || @warn "The southernmost parallel of `target_grid` $φ₁ does not coincide with " *
                          "the closest parallel of the bathymetry grid, $φ₁_data."
    φ₂ ≈ φ₂_data || @warn "The northernmost parallel of `target_grid` $φ₂ does not coincide with " *
                          "the closest parallel of the bathymetry grid, $φ₂_data."

    # Restrict bathymetry _data to region of interest
    λ_data = λ_data[ii]
    φ_data = φ_data[jj]
    h_data = h_data[ii, jj]

    if !isnothing(height_above_water)
        # Overwrite the height of cells above water.
        # This has an impact on reconstruction. Greater height_above_water reduces total
        # wet area by biasing coastal regions to land during bathymetry regridding.
        land = h_data .> 0
        h_data[land] .= height_above_water
    end

    # Build the "native" grid of the bathymetry and make a bathymetry field.
    Nxn = length(λ_data)
    Nyn = length(φ_data)
    Nzn = 1

    native_grid = LatitudeLongitudeGrid(arch;
                                        size = (Nxn, Nyn, Nzn),
                                        latitude = (φ₁, φ₂),
                                        longitude = (λ₁, λ₂),
                                        z = (0, 1),
                                        halo = (10, 10, 1))

    native_h = Field{Center, Center, Nothing}(native_grid)
    set!(native_h, h_data)

    target_h = interpolate_bathymetry_in_passes(native_h, target_grid; 
                                                passes = interpolation_passes,
                                                connected_regions_allowed,
                                                minimum_depth)

    return target_h
end

# Here we can either use `regrid!` (three dimensional version) or `interpolate`
function interpolate_bathymetry_in_passes(native_h, target_grid; 
                                          passes = 10,
                                          connected_regions_allowed = 3,
                                          minimum_depth = 0)
    Nλt, Nφt = Nt = size(target_grid)
    Nλn, Nφn = Nn = size(native_h)
    
    if any(Nt[1:2] .> Nn[1:2]) # We are refining the grid (at least in one direction), more passes will not help!
        target_h = Field{Center, Center, Nothing}(target_grid)
        interpolate!(target_h, native_h)
        return target_h
    end
 
    latitude  = y_domain(target_grid)
    longitude = x_domain(target_grid)

    ΔNλ = floor((Nλn - Nλt) / passes)
    ΔNφ = floor((Nφn - Nφt) / passes)

    Nλ = [Nλn - ΔNλ * pass for pass in 1:passes-1]
    Nφ = [Nφn - ΔNφ * pass for pass in 1:passes-1]

    Nλ = Int[Nλ..., Nλt]
    Nφ = Int[Nφ..., Nφt]

    old_h     = native_h
    TX, TY, _ = topology(target_grid)

    for pass = 1:passes - 1
        new_size = (Nλ[pass], Nφ[pass], 1)

        @debug "pass number $pass with size $new_size"
        new_grid = LatitudeLongitudeGrid(architecture(target_grid),
                                         size = new_size, 
                                     latitude = (latitude[1],  latitude[2]), 
                                    longitude = (longitude[1], longitude[2]), 
                                            z = (0, 1),
                                     topology = (TX, TY, Bounded))

        new_h = Field{Center, Center, Nothing}(new_grid)

        interpolate!(new_h, old_h)
        old_h = new_h
    end

    target_h = Field{Center, Center, Nothing}(target_grid)
    interpolate!(target_h, old_h)

    h_data = Array(interior(target_h, :, :, 1))

    if minimum_depth > 0
        shallow_ocean = h_data .> - minimum_depth
        h_data[shallow_ocean] .= 0
    end

    h_data = remove_lakes!(h_data; connected_regions_allowed)
    set!(target_h, h_data)
    fill_halo_regions!(target_h)

    return target_h
end

"""
    remove_lakes!(h_data; connected_regions_allowed = Inf)

Remove lakes from the bathymetry data stored in `h_data`, by identifying connected regions below sea level 
and removing all but the specified number of largest connected regions (which represent the ocean and 
other possibly disconnected regions like the Mediterranean and the Bering sea).

# Arguments
============

- `h_data`: A 2D array representing the bathymetry data.
- `connected_regions_allowed`: The maximum number of connected regions to keep. 
                               Default is `Inf`, which means all connected regions are kept.

"""
function remove_lakes!(h_data::Field; kw...)
    data = Array(interior(h_data, :, :, 1))
    data = remove_lakes!(data; kw...)

    set!(h_data, data)

    return h_data
end

function remove_lakes!(h_data; connected_regions_allowed = Inf)

    if connected_regions_allowed == Inf
        @info "we are not removing lakes"
        return h_data
    end

    bathtmp = deepcopy(h_data)
    batneg  = zeros(Bool, size(bathtmp)...)
    
    batneg[bathtmp.<0] .= true
    
    connectivity = ImageMorphology.strel(batneg)
    labels = ImageMorphology.label_components(connectivity)
    
    total_elements = zeros(maximum(labels))
    label_elements = zeros(maximum(labels))

    for i in 1:lastindex(total_elements)
        total_elements[i] = length(labels[labels .== i])
        label_elements[i] = i
    end
        
    all_idx = []
    ocean_idx = findfirst(x -> x == maximum(total_elements), total_elements)
    push!(all_idx, label_elements[ocean_idx])
    total_elements = filter((x) -> x != total_elements[ocean_idx], total_elements)
    label_elements = filter((x) -> x != label_elements[ocean_idx], label_elements)

    for _ in 1:connected_regions_allowed
        next_maximum = findfirst(x -> x == maximum(total_elements), total_elements)
        push!(all_idx, label_elements[next_maximum])
        total_elements = filter((x) -> x != total_elements[next_maximum], total_elements)
        label_elements = filter((x) -> x != label_elements[next_maximum], label_elements)
    end
        
    labels = Float64.(labels)

    for i in 1:maximum(labels)
        remove_lake = (&).(Tuple(i != idx for idx in all_idx)...)
        if remove_lake
            labels[labels .== i] .= 1e10 # Fictitious super large number
        end
    end

    # Removing land?
    labels[labels .<  1e10] .= 0 
    labels[labels .== 1e10] .= NaN

    bathtmp .+= labels
    bathtmp[isnan.(bathtmp)] .= 0

    return bathtmp
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

