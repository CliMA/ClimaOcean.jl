include("runtests_setup.jl")

using CUDA
using OrthogonalSphericalShellGrids

@testset "Parameter space test" begin

    for arch in test_architectures
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
                                            major_basins = 1)
        
        grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true)

        free_surface = SplitExplicitFreeSurface(grid; substeps = 20)
        ocean = ocean_simulation(grid; free_surface) 

        backend = JRA55NetCDFBackend(4)
        atmosphere = JRA55_prescribed_atmosphere(arch; backend)
        radiation = Radiation(arch)

        # Fluxes are computed when the model is constructed, so we just test that this works.
        @test begin
            coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
            true
        end
    end
end
