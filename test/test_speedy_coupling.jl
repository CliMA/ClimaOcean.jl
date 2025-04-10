import SpeedyWeather
import ClimaOcean
using Oceananigans
using Dates

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