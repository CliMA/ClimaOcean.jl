
using ClimaOcean
using Oceananigans

# A barotropic ocean
grid = TripolarGrid(size = (700, 400, 1), z = (-10, 0))
bottom_height = regrid_bathymetry(grid)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)
ocean = ocean_simulation(grid)

# A prescribed tidal forcing
Φ = FieldTimeSeries("tidal_forcing.jld2", "Φ")
atmos = PrescribedAtmosphere(Φ.grid, Φ.times; tidal_forcing = Φ)
ocean = ocean_simulation(grid)

barotropic_earth = OceanSeaIceModel(ocean, atmos)
