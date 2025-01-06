using ClimaOcean
using Oceananigans
using SpeedyWeather

arch = GPU()
grid = LatitudeLongitudeGrid(size = (720, 360, 100))