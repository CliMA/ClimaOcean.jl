using Oceananigans.Grids: AbstractGrid
using Oceananigans

import XESMF: Regridder, xesmf_coordinates

const Grids = Union{SpeedyWeather.SpectralGrid, AbstractGrid}

function Regridder(src::Grids, dst::Grids; method::String="bilinear", periodic=true)
    src_coords = xesmf_coordinates(src, Center(), Center(), Center())
    dst_coords = xesmf_coordinates(dst, Center(), Center(), Center())

    return XESMF.Regridder(src_coords, dst_coords; method, periodic)
end

two_dimensionalize(lat::Matrix, lon::Matrix) = lat, lon

function two_dimensionalize(lat::AbstractVector, lon::AbstractVector) 
    Nx  = length(lon)
    Ny  = length(lat)
    lat = repeat(lat', Nx)
    lon = repeat(lon, 1, Ny)
    return lat, lon
end

xesmf_coordinates(grid::SpeedyWeather.SpectralGrid, args...) = xesmf_coordinates(grid.grid, args...)
xesmf_coordinates(grid::SpeedyWeather.RingGrids.AbstractGrid, args...) = 
    throw(ArgumentError("xesmf_coordinates not implemented for grid type $(typeof(grid)), maybe you meant to pass a FullGrid?"))

function xesmf_coordinates(grid::SpeedyWeather.RingGrids.AbstractFullGrid, args...)
    lon  = RingGrids.get_lond(grid)
    lat  = RingGrids.get_latd(grid)
    dlon = lon[2] - lon[1]

    lat_b = [90, 0.5 .* (lat[1:end-1] .+ lat[2:end])..., -90]
    lon_b = [lon[1] - dlon / 2, lon .+ dlon / 2...]

    lat,   lon   = two_dimensionalize(lat,   lon)
    lat_b, lon_b = two_dimensionalize(lat_b, lon_b)

    # Python's xESMF expects 2D arrays with (x, y) coordinates
    # in which y varies in dim=1 and x varies in dim=2
    # therefore we transpose the coordinate matrices
    coords_dictionary = Dict("lat"   => permutedims(lat, (2, 1)),  # φ is latitude
                             "lon"   => permutedims(lon, (2, 1)),  # λ is longitude
                             "lat_b" => permutedims(lat_b, (2, 1)),
                             "lon_b" => permutedims(lon_b, (2, 1)))

    return coords_dictionary
end
