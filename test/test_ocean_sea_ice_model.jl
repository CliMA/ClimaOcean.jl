include("runtests_setup.jl")

using CUDA
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSeaIceModels: above_freezing_ocean_temperature!
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using ClimaSeaIce.SeaIceDynamics
using ClimaSeaIce.Rheologies
using Oceananigans.TimeSteppers: Clock

@inline kernel_melting_temperature(i, j, k, grid, liquidus, S) = @inbounds melting_temperature(liquidus, S[i, j, k])

@testset "Time stepping test" begin
    for dataset in [ECCO4Monthly(), EN4Monthly()]

        start_date = DateTimeProlepticGregorian(1993, 1, 1)
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
                                              interpolation_passes = 5,
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

            sea_ice  = sea_ice_simulation(grid, ocean; advection=WENO(order=7))
            liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus

            # Set the ocean temperature and salinity
            set!(ocean.model, T=temperature_metadata[1], S=salinity_metadata[1])
            above_freezing_ocean_temperature!(ocean, grid, sea_ice)

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

@testset "Prescribed atmosphere with DateTime clock" begin
    start_time = Dates.DateTime(1993, 2, 1)
    Δt_seconds = 3050.0
    Δt_period = Dates.Second(3050)
    steps = 5
    stop_time = start_time + steps * Δt_period
    times = [start_time + m * Δt_period for m in 0:steps]

    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size = (1, 1, 20),
                                     longitude = (10, 11),
                                     latitude = (20.2, 21.2),
                                     z = (-200, 0),
                                     topology = (Bounded, Bounded, Bounded))

        ocean_clock = Clock(time=start_time)
        ocean = ocean_simulation(grid; clock=ocean_clock, Δt=Δt_seconds, tracers=(:T, :S), verbose=false, closure=nothing)
        atmosphere_clock = Clock(time=start_time)
        atmosphere = PrescribedAtmosphere(grid, times; clock=atmosphere_clock)
        radiation = Radiation(arch)
        coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

        for _ in 1:steps
            time_step!(coupled_model, Δt_seconds)
        end

        @test coupled_model.clock.time == stop_time
        @test coupled_model.clock.iteration == steps
        @test atmosphere.clock.time == coupled_model.clock.time
        @test coupled_model.clock.time isa Dates.DateTime

        ocean_clock = Clock(time=start_time)
        ocean = ocean_simulation(grid; clock=ocean_clock, Δt=Δt_seconds, tracers=(:T, :S), verbose=false, closure=nothing)
        atmosphere_clock = Clock(time=start_time)
        atmosphere = PrescribedAtmosphere(grid, times; clock=atmosphere_clock)
        radiation = Radiation(arch)
        coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)
        coupled_simulation = Simulation(coupled_model; Δt=Δt_seconds, stop_time)

        run!(coupled_simulation)

        @test coupled_simulation.model.clock.time == stop_time
        @test coupled_simulation.model.clock.iteration == steps
        @test coupled_simulation.model.atmosphere.clock.time == coupled_model.clock.time
        @test coupled_simulation.model.clock.time isa Dates.DateTime
    end
end
