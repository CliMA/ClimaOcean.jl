include("runtests_setup.jl")

using Glob
using Oceananigans.OutputWriters: Checkpointer

@testset "OceanSeaIceModel checkpointing" begin
    for arch in test_architectures
        arch = CPU()
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
            set!(ocean.model, T=20, S=35, u=0.01)
            sea_ice = sea_ice_simulation(grid, ocean)
            set!(sea_ice.model, h=hi, ℵ=hi)

            # Create atmosphere and radiation
            backend = JRA55NetCDFBackend(4)
            atmosphere = JRA55PrescribedAtmosphere(arch; backend)

            return OceanSeaIceModel(ocean, sea_ice; atmosphere)
        end

        # Reference run: 6 iterations continuously
        model = make_coupled_model(grid)
        simulation = Simulation(model, Δt=60, stop_iteration=6)
        run!(simulation)

        # Store reference states
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
        simulation.output_writers[:checkpointer] = Checkpointer(model;
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
        run!(simulation)

        # Compare final states
        T  = Array(interior(model.ocean.model.tracers.T))
        S  = Array(interior(model.ocean.model.tracers.S))
        u  = Array(interior(model.ocean.model.velocities.u))
        v  = Array(interior(model.ocean.model.velocities.v))
        h  = Array(interior(model.sea_ice.model.ice_thickness))
        ui = Array(interior(model.sea_ice.model.velocities.u))
        vi = Array(interior(model.sea_ice.model.velocities.v))

        @test all(T  .≈ ref_T)
        @test all(S  .≈ ref_S)
        @test all(u  .≈ ref_u)
        @test all(v  .≈ ref_v)
        @test all(h  .≈ ref_h)
        @test all(ui .≈ ref_ui)
        @test all(vi .≈ ref_vi)
        @test model.clock.time ≈ ref_time
        @test model.clock.iteration == ref_iteration

        # Cleanup
        rm.(glob("$(prefix)_iteration*.jld2"), force=true)
    end
end
