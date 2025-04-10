import SpeedyWeather
import ClimaOcean
using Oceananigans
using Dates
using Test

ClimaOceanSpeedyWeatherExt = Base.get_extension(ClimaOcean, :ClimaOceanSpeedyWeatherExt)
@test !isnothing(ClimaOceanSpeedyWeatherExt)

spectral_grid = SpeedyWeather.SpectralGrid(trunc=31, nlayers=10)
oceananigans_grid = ClimaOceanSpeedyWeatherExt.clenshaw_latitude_longitude_grid(CPU(), spectral_grid)
λd_sw, φd_sw = SpeedyWeather.RingGrids.get_londlatds(spectral_grid.Grid, spectral_grid.nlat_half)
λd_oc = λnodes(oceananigans_grid, Center(), Center(), Center())
φd_oc = φnodes(oceananigans_grid, Center(), Center(), Center())

@test λd_sw ≈ λd_oc
@test φd_sw ≈ φd_oc

# SpeedyWeather.RingGrids.update_locator!(interpolator, londs, latds)

#=

arch = CPU()
Nx = 180
Ny = 85
Nz = 10
longitude = (0, 360)
latitude = (-85, 85)
halo = (7, 7, 7)
ocean_grid = LatitudeLongitudeGrid(arch; size=(Nx, Ny, Nz), halo, longitude, latitude, z=(-100, 0))
ocean = ClimaOcean.ocean_simulation(ocean_grid)

atmos_grid = SpeedyWeather.SpectralGrid(trunc=31, nlayers=10)
atmosphere = ClimaOcean.atmosphere_simulation(atmos_grid)
SpeedyWeather.set!(atmosphere.model.time_stepping, Δt=Minute(10))

# EarthSystemModel
# const OceanSeaIceModel = 
# const SeaIceOnlyModel = 
# const LandOnlyModel = 

model = ClimaOcean.OceanSeaIceModel(ocean; atmosphere)
=#