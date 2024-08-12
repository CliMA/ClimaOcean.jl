using Oceananigans
include("runtests_setup.jl")

using OrthogonalSphericalShellGrids

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

free_surface = SplitExplicitFreeSurface(grid; substeps = 20)
ocean = ocean_simulation(grid; free_surface) 
model = ocean.model

backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation(arch)

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

# Fluxes are computed when the model is constructed, so we just test that this works.
@test begin
    coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
    true
end

