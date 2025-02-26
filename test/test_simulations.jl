include("runtests_setup.jl")

using CUDA
using OrthogonalSphericalShellGrids

@testset "GPU time stepping test" begin

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
        atmosphere = JRA55PrescribedAtmosphere(arch; backend)
        radiation = Radiation(arch)

        # Fluxes are computed when the model is constructed, so we just test that this works.
        @test begin
            coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
            time_step!(coupled_model, 1)
            true
        end
    end
end
