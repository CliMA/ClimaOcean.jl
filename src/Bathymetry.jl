module Bathymetry

using ..DataWrangling: download_progress

using Oceananigans
using Oceananigans.Architectures: architecture
using Oceananigans.Grids: halo_size, λnodes
using Oceananigans.Utils: pretty_filesize

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

Arguments
=========

Keyword arguments
=================
"""
function regrid_bathymetry(target_grid;
                           height_above_water = nothing,
                           minimum_depth = 0,
                           dir = joinpath(@__DIR__, "..", "data"),
                           url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf", 
                           filename = "ETOPO_2022_v1_60s_N90W180_surface.nc")

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

    φ_data = dataset["lat"][:]
    λ_data = dataset["lon"][:]
    h_data = dataset["z"][:, :]

    close(dataset)

    # Diagnose target grid information
    Nxt, Nyt, Nzt = size(target_grid)
    arch = architecture(target_grid)
    λt, φt, zt = nodes(target_grid, Face(), Face(), Face())

    λ₁ = λt[1]
    λ₂ = λt[Nxt]

    φ₁ = φt[1]
    φ₂ = φt[Nyt]

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

    if minimum_depth > 0
        ocean = h_data .<= 0
        h_data[ocean] .= min.(-minimum_depth, h_data[ocean])
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
                                        halo = halo_size(target_grid))

    native_h = Field{Center, Center, Nothing}(native_grid)
    set!(native_h, h_data)

    Nxi = Nxt
    Nyi = Nyn
    Nzi = Nzn

    if parent(parent(λt)) isa StepRangeLen # longitude is equispaced
        longitude = (λ₁, λ₂)
    else
        longitude = λnodes(target_grid, Face(), Center(), Center())
    end

    intermediate_grid = LatitudeLongitudeGrid(arch;
                                              size = (Nxi, Nyi, Nzi),
                                              latitude = (φ₁, φ₂),
                                              longitude,
                                              z = (0, 1),
                                              halo = halo_size(target_grid))

    intermediate_h =  Field{Center, Center, Nothing}(intermediate_grid)
    regrid!(intermediate_h, native_h)

    one_degree_h = Field{Center, Center, Nothing}(target_grid)
    regrid!(one_degree_h, intermediate_h)

    return one_degree_h
end

end # module

