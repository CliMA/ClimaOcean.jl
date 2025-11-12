using SpeedyWeather, XESMF
using ClimaOcean
using Oceananigans
using Dates
using Test

ClimaOceanSpeedyWeatherExt = Base.get_extension(ClimaOcean, :ClimaOceanSpeedyWeatherExt)
@test !isnothing(ClimaOceanSpeedyWeatherExt)

spectral_grid = SpeedyWeather.SpectralGrid(trunc=51, nlayers=3, Grid=FullClenshawGrid)
oceananigans_grid = LatitudeLongitudeGrid(Oceananigans.CPU(); size=(200, 100, 1), latitude=(-80, 80), longitude=(0, 360), z = (0, 1))

ocean = ClimaOcean.OceanSimulations.ocean_simulation(oceananigans_grid; momentum_advection=nothing, tracer_advection=nothing, closure=nothing)
Oceananigans.set!(ocean.model, T=EN4Metadatum(:temperature), S=EN4Metadatum(:salinity))

atmos = ClimaOcean.atmosphere_simulation(spectral_grid)

radiation   = Radiation(ocean_emissivity=0.0, sea_ice_emissivity=0.0)
earth_model = OceanSeaIceModel(ocean; atmosphere=atmos, radiation)

Qca = atmos.prognostic_variables.ocean.sensible_heat_flux.data
Mva = atmos.prognostic_variables.ocean.surface_humidity_flux.data

@test !(all(Qca .== 0.0))
@test !(all(Mva .== 0.0))