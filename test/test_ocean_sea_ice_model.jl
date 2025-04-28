include("runtests_setup.jl")

using CUDA
using Oceananigans.Units
using Oceananigans: prognostic_fields
using Oceananigans.TimeSteppers: QuasiAdamsBashforth2TimeStepper
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSeaIceModels: above_freezing_ocean_temperature!, OceanSeaIceModel
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies

@inline kernel_melting_temperature(i, j, k, grid, liquidus, S) = @inbounds melting_temperature(liquidus, S[i, j, k])

@testset "Time stepping test" begin
    for dataset in test_datasets

        start = DateTimeProlepticGregorian(1993, 1, 1)
        time_resolution = dataset isa ECCO2Daily ? Day(1) : Month(1)
        end_date = DateTimeProlepticGregorian(1993, 2, 1)
        dates = start_date : time_resolution : end_date

        temperature_metadata = Metadata(:temperature; dataset, dates)
        salinity_metadata    = Metadata(:salinity; dataset, dates)

        for arch in test_architectures

            A = typeof(arch)

            @info "Testing timestepping with $(typeof(dataset)) on $A"

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
            @test all(T .≥ Tm)

            # Fluxes are computed when the model is constructed, so we just test that this works.
            # And that we can time step with sea ice
            @test begin
                coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
                time_step!(coupled_model, 1)
                true
            end
        end
    end
end

"""
    test_ocean_sea_ice_model_checkpointer(grid)

Return two ocean-sea ice coupled models to be used for testing purposes.
"""
function test_ocean_sea_ice_model_checkpointer(grid)
    ocean = ocean_simulation(grid, closure=nothing, free_surface = ExplicitFreeSurface())
    atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(4))

    true_coupled_model = OceanSeaIceModel(ocean; atmosphere)
    # true_coupled_model = OceanSeaIceModel(ocean; atmosphere=nothing, radiation=nothing)

    test_coupled_model = deepcopy(true_coupled_model)

    return true_coupled_model, test_coupled_model
end

function test_model_equality(test_model::OceanSeaIceModel, true_model::OceanSeaIceModel)
    test_clock_equality(test_model.clock, true_model.clock)

    # test equality for components
    test_model_equality(test_model.atmosphere, true_model.atmosphere)
    test_model_equality(test_model.ocean.model, true_model.ocean.model)

    return nothing
end

test_model_equality(test_model::PrescribedAtmosphere, true_model::PrescribedAtmosphere) =
    test_clock_equality(test_model.clock, true_model.clock)

function test_clock_equality(clock1::Clock, clock2::Clock)
    @test clock1.iteration == clock2.iteration
    @test clock1.time == clock2.time
    @test clock1.last_Δt == clock2.last_Δt
    return nothing
end

function test_model_equality(test_model::HydrostaticFreeSurfaceModel, true_model::HydrostaticFreeSurfaceModel)

    test_clock_equality(test_model.clock, true_model.clock)

    CUDA.@allowscalar begin
        test_model_fields = prognostic_fields(test_model)
        true_model_fields = prognostic_fields(true_model)
        field_names = keys(test_model_fields)

        for name in field_names
            @show name
            @show all(test_model_fields[name].data .≈ true_model_fields[name].data)

            if test_model.timestepper isa QuasiAdamsBashforth2TimeStepper
                if name ∈ keys(test_model.timestepper.Gⁿ)
                    @show name
                    # @show test_model.timestepper.Gⁿ[name]
                    # @show true_model.timestepper.Gⁿ[name]
                    @show all(test_model.timestepper.Gⁿ[name].data .≈ true_model.timestepper.Gⁿ[name].data)
                    @show all(test_model.timestepper.G⁻[name].data .≈ true_model.timestepper.G⁻[name].data)
                end
            end
        end
    end

    return nothing
end

function run_checkpointer_tests(true_model, test_model, Δt)

    true_simulation = Simulation(true_model; Δt, stop_iteration = 4)
    true_simulation.output_writers[:checkpoint] =
        Checkpointer(true_model, schedule=IterationInterval(4), overwrite_existing=true)

    run!(true_simulation)

    checkpointed_model = deepcopy(true_simulation.model)

    true_simulation.stop_iteration = 7
    run!(true_simulation)

    #####
    ##### Test `set!(model, checkpoint_file)`
    #####

    @info "Testing `set!(model, checkpoint_file)`"

    set!(test_model, "checkpoint_iteration4.jld2")

    test_model_equality(test_model, checkpointed_model)

    #=
    #####
    ##### Test pickup from explicit checkpoint path
    #####

    @info "Testing pickup from explicit checkpoint path"

    test_simulation = Simulation(test_model; Δt, stop_iteration=7)
    run!(test_simulation, pickup="checkpoint_iteration0.jld2")

    test_model_equality(test_model, true_model)

    run!(test_simulation, pickup="checkpoint_iteration4.jld2")

    @info "Testing model equality when running with pickup=checkpoint_iteration4.jld2."

    test_model_equality(test_model, true_model)

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation.model, true_simulation.model)

    #####
    ##### Test `run!(sim, pickup=true)
    #####

    @info "Testing run!(sim, pickup=true)"

    # Pickup using existing checkpointer
    test_simulation.output_writers[:checkpoint] =
        Checkpointer(test_model, schedule=IterationInterval(4), overwrite_existing=true)

    run!(test_simulation, pickup=true)
    @info "    Testing model equality when running with pickup=true."

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation.model, true_simulation.model)

    run!(test_simulation, pickup=0)
    @info "    Testing model equality when running with pickup=0."

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation.model, true_simulation.model)

    run!(test_simulation, pickup=4)
    @info "    Testing model equality when running with pickup=4."

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation.model, true_simulation.model)

    =#

    for file in Tuple("checkpoint_iteration$i.jld2" for i in 0:4:10)
        rm(file, force=true)
    end
end

#=
@testset "Checkpointer tests" begin
    for arch in test_architectures

        Δt = 1.2hours

        Nx, Ny, Nz = 40, 25, 10

        grid = LatitudeLongitudeGrid(arch;
                                     size = (Nx, Ny, Nz),
                                     halo = (7, 7, 7),
                                     z = (-6000, 0),
                                     latitude  = (-75, 75),
                                     longitude = (0, 360))

        true_coupled_model, test_coupled_model = test_ocean_sea_ice_model_checkpointer(grid)

        run_checkpointer_tests(true_coupled_model, test_coupled_model, Δt)

        #=
        true_simulation = testbed_coupled_simulation(true_coupled_model; Δt, checkpoint_iteration_interval,
                                                     stop_iteration=intermediate_iteration)

        test_simulation = testbed_coupled_simulation(test_coupled_model; Δt, checkpoint_iteration_interval
                                                     stop_iteration = final_iteration)

        run!(true_simulation) # up to intermediate_iteration

        checkpointed_model = deepcopy(true_simulation.model)

        set!(test_model, "checkpoint_iteration5.jld2")

        @test test_model.clock.iteration == checkpointed_model.clock.iteration
        @test test_model.clock.time == checkpointed_model.clock.time

        run!(new_simulation, pickup=true)

        # ensure all clocks all at same time and iteration
        for clock in (simulation.model.clock,
                      simulation.model.atmosphere.clock,
                      simulation.model.ocean.model.clock)

            @test clock.iteration ≈ intermediate_iteration
            @test clock.time ≈ intermediate_iteration * Δt
        end

        for clock in (new_simulation.model.clock,
                      new_simulation.model.atmosphere.clock,
                      new_simulation.model.ocean.model.clock)

            @test clock.iteration ≈ final_iteration
            @test clock.time ≈ final_iteration * Δt
        end
        =#
    end

end
=#
