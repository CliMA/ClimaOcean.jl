include("runtests_setup.jl")

using CUDA
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSeaIceModels: above_freezing_ocean_temperature!
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies

@inline kernel_melting_temperature(i, j, k, grid, liquidus, S) = @inbounds melting_temperature(liquidus, S[i, j, k])

@testset "Time stepping test" begin

    for arch in test_architectures

        λ★, φ★ = 35.1, 50.1

        grid = RectilinearGrid(arch, size = 200, x = λ★, y = φ★,
                               z = (-400, 0), topology = (Flat, Flat, Bounded))

        ocean = ocean_simulation(grid)
        data = Int[]
        pushdata(sim) = push!(data, iteration(sim))
        add_callback!(ocean, pushdata)
        backend = JRA55NetCDFBackend(4)
        atmosphere = JRA55PrescribedAtmosphere(arch; backend)
        radiation = Radiation(arch)
        coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
        Δt = 60
        for n = 1:3
            time_step!(coupled_model, Δt)
        end
        @test data == [0, 1, 2, 3]

        # TODO: do the same for a SeaIceSimulation, and eventually prognostic Atmos

        #####
        ##### Ocean and prescribed atmosphere
        #####

        grid = TripolarGrid(arch;
                            size = (50, 50, 10),
                            halo = (7, 7, 7),
                            z = (-6000, 0))

        bottom_height = regrid_bathymetry(grid; 
                                          minimum_depth = 10,
                                          interpolation_passes = 20,
                                          major_basins = 1)
        
        grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

        free_surface = SplitExplicitFreeSurface(grid; substeps=20)
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

        # Adding a sea ice model to the coupled model
        τua = Field{Face, Center, Nothing}(grid) 
        τva = Field{Center, Face, Nothing}(grid)

        dynamics = SeaIceMomentumEquation(grid; 
                                          coriolis = ocean.model.coriolis,
                                          top_momentum_stress = (u=τua, v=τva),
                                          rheology = ElastoViscoPlasticRheology(),
                                          solver = SplitExplicitSolver(120))

        sea_ice  = sea_ice_simulation(grid; dynamics, advection=WENO(order=7)) 
        liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus
        
        # Set the ocean temperature and salinity
        set!(ocean.model, T=temperature_metadata[1], S=salinity_metadata[1])
        above_freezing_ocean_temperature!(ocean, sea_ice)

        # Test that ocean temperatures are above freezing
        T = on_architecture(CPU(), ocean.model.tracers.T)
        S = on_architecture(CPU(), ocean.model.tracers.S)

        Tm = KernelFunctionOperation{Center, Center, Center}(kernel_melting_temperature, grid, liquidus, S)
        @test all(T .>= Tm)

        # Fluxes are computed when the model is constructed, so we just test that this works.
        # And that we can time step with sea ice
        @test begin
            coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
            time_step!(coupled_model, 1)
            true
        end
    end

    @testset "Test a prescribed ocean model" begin 
        grid = LatitudeLongitudeGrid(size = (10, 10, 10), latitude = (20, 30), longitude = (40, 50), z = (-100, 0))
        T = ECCOFieldTimeSeries(:temperature, grid; time_indices_in_memory = 2)
        S = ECCOFieldTimeSeries(:salinity, grid; time_indices_in_memory = 2)
        ocean_model = PrescribedOcean((; T, S); grid)

        time_step!(ocean_model, T.times[2])

        @test ocean_model.clock.time == T.times[2]
        @test ocean_model.tracers.T.data == T[2].data
        @test ocean_model.tracers.S.data == S[2].data

        time_step!(ocean_model, 10days)

        @test ocean_model.clock.time == 10days
        @test T.backend.start == 2
        @test ocean_model.tracers.T.data != T[3].data 
    end
end

