module ClimaOceanSpeedyWeatherExt

using OffsetArrays
using KernelAbstractions
using Statistics

import SpeedyWeather 
import ClimaOcean 
import Oceananigans 
import SpeedyWeather.RingGrids

include("speedy_atmosphere_simulations.jl")

# TODO: Remove this when we have a proper interface
# for regridding between tripolar and latitude-longitude grids
include("bilinear_interpolator.jl")
include("speedy_weather_exchanger.jl")

end # module ClimaOceanSpeedyWeatherExt