using Oceananigans
using ClimaOcean
using OrthogonalSphericalShellGrids

start_time = time_ns()
arch = GPU() 
grid = TripolarGrid(arch; 
                    size = (50, 50, 10), 
                    halo = (7, 7, 7), 
                    z = (-6000, 0), 
                    first_pole_longitude = 75,
                    north_poles_latitude = 55)

bottom_height = retrieve_bathymetry(grid; 
                                    minimum_depth = 10,
                                    dir = "./",
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 0)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true)                           

elapsed = 1e-9 * (time_ns() - start_time)
@info "Grid / bathymetry construction time: " * prettytime(elapsed)

start_time = time_ns()
free_surface = SplitExplicitFreeSurface(grid; substeps = 20)
ocean = ocean_simulation(grid; free_surface) 
model = ocean.model
@info "Ocean simulation construction time: " * prettytime(elapsed)

start_time = time_ns()
backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation(arch)

elapsed = 1e-9 * (time_ns() - start_time)
@info "Atmosphere construction time: " * prettytime(elapsed)

# Fluxes are computed when the model is constructed, so we just test that this works.
start_time = time_ns()
sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

elapsed = 1e-9 * (time_ns() - start_time)
@info "Coupled model construction time: " * prettytime(elapsed)

start_time = time_ns()
time_step!(coupled_model, 1)
elapsed = 1e-9 * (time_ns() - start_time)
@info "One time step time: " * prettytime(elapsed)
