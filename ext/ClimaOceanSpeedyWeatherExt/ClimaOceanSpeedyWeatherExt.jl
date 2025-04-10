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

# include("conservative_regridding.jl")
include("speedy_atmosphere_simulations.jl")
include("speedy_weather_fluxes.jl")

end # module ClimaOceanSpeedyWeatherExt