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
end

"""
    testbed_coupled_simulation(grid; stop_iteration=8)

Return a test-bed coupled simulation with a Checkpointer.
"""
function testbed_coupled_simulation(grid; stop_iteration=8)
    ocean = ocean_simulation(grid)

    atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(4))

    coupled_model = OceanSeaIceModel(ocean; atmosphere)

    simulation = Simulation(coupled_model; Δt=10, stop_iteration)

    simulation.output_writers[:checkpoint] = Checkpointer(ocean.model;
                                                          schedule = IterationInterval(3),
                                                          prefix = "checkpointer",
                                                          dir = ".",
                                                          verbose = true,
                                                          overwrite_existing = true)

    return simulation
end

@testset "Checkpointer test" begin
    for arch in test_architectures

        Nx, Ny, Nz = 40, 25, 10

        grid = LatitudeLongitudeGrid(arch;
                                     size = (Nx, Ny, Nz),
                                     halo = (7, 7, 7),
                                     z = (-6000, 0),
                                     latitude  = (-75, 75),
                                     longitude = (0, 360))

        simulation = testbed_coupled_simulation(grid; stop_iteration=8)

        run!(simulation)

        # create a new simulation and pick up from the last checkpointer
        new_simulation = testbed_coupled_simulation(grid; stop_iteration=13)

        run!(new_simulation, pickup=true)

        # ensure the ocean, atmosphere, and coupled model are all at same time and iteration
        clock = new_simulation.model.ocean.model.clock
        @test new_simulation.model.atmosphere.clock.iteration ≈ clock.iteration
        @test new_simulation.model.atmosphere.clock.time ≈ clock.time
        @test new_simulation.model.clock.iteration ≈ clock.iteration
        @test new_simulation.model.clock.time ≈ clock.time
    end
end
