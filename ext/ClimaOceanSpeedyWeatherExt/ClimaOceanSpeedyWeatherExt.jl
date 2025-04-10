module ClimaOceanSpeedyWeatherExt

using OffsetArrays
using KernelAbstractions
using Statistics

import SpeedyWeather 
import ClimaOcean 
import Oceananigans 

include("conservative_regridding.jl")
include("speedy_weather_fluxes.jl")

end # module ClimaOceanSpeedyWeatherExt