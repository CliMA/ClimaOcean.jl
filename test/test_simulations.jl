include("runtests_setup.jl")

using CUDA
using OrthogonalSphericalShellGrids
using ClimaOcean.OceanSeaIceModels: adjust_freezing_ocean_temperature!
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

@testset "GPU time stepping test" begin

    for arch in test_architectures

        #####
        ##### Ocean and prescribed atmosphere
        #####

        grid = TripolarGrid(arch;
                            size = (50, 50, 10),
                            halo = (7, 7, 7),
                            z = (-6000, 0)

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

        #####
        ##### Coupled ocean-sea ice and prescribed atmosphere
        #####

        sea_ice_grid = TripolarGrid(arch; size=(50, 50, 1), z = (-10, 0))
        sea_ice_grid = ImmersedBoundaryGrid(sea_ice_grid, GridFittedBottom(bottom_height))

        sea_ice  = sea_ice_simulation(sea_ice_grid) 
        liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus
        
        # Set the ocean temperature and salinity
        set!(ocean.model, T=temperature_metadata[1], S=salinity_metadata[1])

        adjust_freezing_ocean_temperature!(ocean, sea_ice)

        # Test that ocean temperatures are above freezing
        T = on_architecture(CPU(), ocean.model.T)
        S = on_architecture(CPU(), ocean.model.S)

        @inline pointwise_melting_T(i, j, k, grid, liquidus, S) = @inbounds melting_temperature(liquidus, S[i, j, k])

        Tm = KernelFunctionOperation{Center, Center, Center}(pointwise_melting_T, grid, S)

        @test all(T .> Tm)

        # Fluxes are computed when the model is constructed, so we just test that this works.
        # And that we can time step with sea ice
        @test begin
            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
            true
        end
    end
end
