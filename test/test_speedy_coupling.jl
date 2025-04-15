using SpeedyWeather
using ClimaOcean
using Oceananigans
using Dates
using Test

ClimaOceanSpeedyWeatherExt = Base.get_extension(ClimaOcean, :ClimaOceanSpeedyWeatherExt)
@test !isnothing(ClimaOceanSpeedyWeatherExt)

spectral_grid = SpeedyWeather.SpectralGrid(trunc=31, nlayers=10)
oceananigans_grid = LatitudeLongitudeGrid(Oceananigans.CPU(); size=(200, 100, 1), latitude=(-80, 80), longitude=(0, 360), z = (0, 1))

ocean = ClimaOcean.OceanSimulations.ocean_simulation(oceananigans_grid; momentum_advection=nothing, tracer_advection=nothing, closure=nothing)
atmos = ClimaOcean.atmosphere_simulation(spectral_grid)
earth = OceanSeaIceModel(ocean; atmosphere=atmos)

fluxes = earth.interfaces.atmosphere_ocean_interface.fluxes

Oceananigans.set!(fluxes.sensible_heat, (x, y) -> exp(- y^2 / 10^2))

# Regridding to speedy
ClimaOcean.OceanSeaIceModel.InterfaceComputations.compute_net_atmosphere_fluxes!(earth)