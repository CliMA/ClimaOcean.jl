module ClimaOceanSpeedyWeatherExt

using OffsetArrays
using KernelAbstractions
using Statistics

import SpeedyWeather 
import ClimaOcean 
import Oceananigans 

function clenshaw_latitude_longitude_grid(arch, spectral_grid)
    grid = LatitudeLongitudeGrid(arch;
                                 size = (360, 179, 1),
                                 latitude = (-89.5, 89.5),
                                 longitude = (0, 360),
                                 z = (0, 1))
    return grid
end

# Put it here for the moment, 
# but use the one from Speecy before merging
function get_faces(
    Grid::Type{<:SpeedyWeather.AbstractGridArray},
    nlat_half::Integer;
    add_nan::Bool = false,
)
    npoints = SpeedyWeather.RingGrids.get_npoints2D(Grid, nlat_half)

    # vertex east, south, west, north (i.e. clockwise for every grid point)
    E, S, W, N = SpeedyWeather.RingGrids.get_vertices(Grid, nlat_half)

    @boundscheck size(N) == size(W) == size(S) == size(E) || throw(BoundsError("Vertices must have the same size"))
    @boundscheck size(N) == (2, npoints) || throw(BoundsError("Number of vertices and npoints do not agree"))

    # number of vertices = 4, 5 to close the polygon, 6 to add a nan
    # to prevent grid lines to be drawn between cells
    nvertices = add_nan ? 6 : 5

    # allocate faces as Point2{Float64} so that no data copy has to be made in Makie
    faces = Matrix{NTuple{2, Float64}}(undef, nvertices, npoints)

    @inbounds for ij in 1:npoints
        faces[1, ij] = (E[1, ij], E[2, ij])  # clockwise
        faces[2, ij] = (S[1, ij], S[2, ij])
        faces[3, ij] = (W[1, ij], W[2, ij])
        faces[4, ij] = (N[1, ij], N[2, ij])
        faces[5, ij] = (E[1, ij], E[2, ij])  # back to east to close the polygon        
    end

    if add_nan  # add a NaN to separate grid cells
        for ij in 1:npoints
            faces[6, ij] = (NaN, NaN)
        end
    end

    return faces
end

get_faces(grid::SpeedyWeather.AbstractGridArray; kwargs...) = get_faces(typeof(grid), grid.nlat_half; kwargs...)

include("speedy_atmosphere_simulations.jl")
include("speedy_weather_exchanger.jl")
# include("speedy_weather_fluxes.jl")

end # module ClimaOceanSpeedyWeatherExt