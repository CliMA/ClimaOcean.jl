include("runtests_setup.jl")

using Glob
using Oceananigans.OutputWriters: Checkpointer

@testset "OceanSeaIceModel checkpointing" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing OceanSeaIceModel checkpointing on $A"

        # Create a minimal grid
        grid = LatitudeLongitudeGrid(arch;
                                     size = (20, 20, 4),
                                     z = (-100, 0),
                                     latitude = (-80, 80),
                                     longitude = (0, 360),
                                     halo = (6, 6, 6))

        # Helper function to create fresh simulations
        function make_coupled_model(grid)
            @inline hi(λ, φ) = φ > 70 || φ < -70

            # Create ocean and Sea ice
            ocean = ocean_simulation(grid, closure=nothing)
            set!(ocean.model, T=20, S=35, u=0.01, v=-0.005)
            sea_ice = sea_ice_simulation(grid, ocean)
            set!(sea_ice.model, h=hi, ℵ=hi)

            # Create atmosphere and radiation
            backend = JRA55NetCDFBackend(4)
            atmosphere = JRA55PrescribedAtmosphere(arch; backend)

            return OceanSeaIceModel(ocean, sea_ice; atmosphere)
        end

        # Reference run: 3 iterations, then continue to 6
        # (This matches the checkpointed workflow where we create a new Simulation
        # after iteration 3, which is what happens during checkpoint restore)
        model = make_coupled_model(grid)
        simulation = Simulation(model, Δt=60, stop_iteration=3)
        run!(simulation)

        # Continue on the same model (simulates what happens after checkpoint restore)
        simulation = Simulation(model, Δt=60, stop_iteration=6)
        run!(simulation)

        # Store reference states at iteration 6
        ref_T  = Array(interior(model.ocean.model.tracers.T))
        ref_S  = Array(interior(model.ocean.model.tracers.S))
        ref_u  = Array(interior(model.ocean.model.velocities.u))
        ref_v  = Array(interior(model.ocean.model.velocities.v))
        ref_h  = Array(interior(model.sea_ice.model.ice_thickness))
        ref_ui = Array(interior(model.sea_ice.model.velocities.u))
        ref_vi = Array(interior(model.sea_ice.model.velocities.v))
        ref_time = model.clock.time
        ref_iteration = model.clock.iteration

        # Checkpointed run: 3 iterations, then checkpoint
        model = make_coupled_model(grid)
        simulation = Simulation(model, Δt=60, stop_iteration=3)

        prefix = "osm_checkpointer_test_$(typeof(arch))"
        simulation.output_writers[:checkpointer] = Checkpointer(simulation.model;
                                                                schedule = IterationInterval(3),
                                                                prefix = prefix)

        run!(simulation)

        @test isfile("$(prefix)_iteration3.jld2")

        # Create new model and restore from checkpoint
        model = make_coupled_model(grid)
        simulation = Simulation(model, Δt=60, stop_iteration=6)

        simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                                schedule = IterationInterval(3),
                                                                prefix = prefix)

        set!(simulation; checkpoint=:latest)

        @test simulation.model.clock.iteration == 3

        set!(simulation; iteration=3)

        @test simulation.model.clock.iteration == 3

        run!(simulation)

        # Compare final states
        T  = Array(interior(model.ocean.model.tracers.T))
        S  = Array(interior(model.ocean.model.tracers.S))
        u  = Array(interior(model.ocean.model.velocities.u))
        v  = Array(interior(model.ocean.model.velocities.v))
        h  = Array(interior(model.sea_ice.model.ice_thickness))
        ui = Array(interior(model.sea_ice.model.velocities.u))
        vi = Array(interior(model.sea_ice.model.velocities.v))

        # Checkpoint restore produces results within numerical precision of the iterative solvers.
        # The split-explicit (ocean) and EVP (sea ice) solvers accumulate floating point
        # differences during substepping, even with identical initial conditions.
        @test T ≈ ref_T rtol=1e-13
        @test S ≈ ref_S rtol=1e-13
        @test h ≈ ref_h rtol=1e-13
        @test u ≈ ref_u rtol=1e-10  # split-explicit solver precision
        @test v ≈ ref_v rtol=1e-10  # split-explicit solver precision
        @test ui ≈ ref_ui rtol=1e-10  # EVP solver precision
        @test vi ≈ ref_vi rtol=1e-10  # EVP solver precision
        @test model.clock.time == ref_time
        @test model.clock.iteration == ref_iteration

        # Cleanup
        rm.(glob("$(prefix)_iteration*.jld2"), force=true)
    end
end
