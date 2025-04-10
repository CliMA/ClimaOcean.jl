import SpeedyWeather
import ClimaOcean
using Oceananigans

arch = CPU()
Nx = 180
Ny = 85
Nz = 10
longitude = (0, 360)
latitude = (-85, 85)
grid = LatitudeLongitudeGrid(arch; size=(Nx, Ny, Nz), longitude, latitude, z=(-100, 0))
ocean = ClimaOcean.ocean_simulation(grid)