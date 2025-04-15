import SpeedyWeather
import ClimaOcean
using Oceananigans
using Dates
using Test

ClimaOceanSpeedyWeatherExt = Base.get_extension(ClimaOcean, :ClimaOceanSpeedyWeatherExt)
@test !isnothing(ClimaOceanSpeedyWeatherExt)

spectral_grid = SpeedyWeather.SpectralGrid(trunc=31, nlayers=10)
oceananigans_grid = LatitudeLongitudeGrid(CPU(); size=(200, 100, 1), latitude=(-80, 80), longitude=(0, 360), z = (0, 1))

ocean = ocean_simulation(oceananigans_grid; momentum_advection=nothing, tracer_advection=nothing, closure=nothing)
atmos = atmosphere_simulation(spectral_grid)
earth = OceanSeaIceModel(ocean; atmosphere=atmos)