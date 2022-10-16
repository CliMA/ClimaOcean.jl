using Oceananigans

using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation

simulation = one_degree_near_global_simulation()

start_time = [time_ns()]

function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    u = sim.model.velocities.u
    w = sim.model.velocities.w

    intw  = Array(interior(w))
    max_w = findmax(intw)

    mw = max_w[1]
    iw = max_w[2]

    @info @sprintf("Time: % 12s, iteration: %d, max(|u|): %.2e ms⁻¹, wmax: %.2e , loc: (%d, %d, %d), wall time: %s", 
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration, maximum(abs, u), mw, iw[1], iw[2], iw[3], 
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S

simulation.output_writers[:fields] =
    JLD2OutputWriter(model, (; u, v, w, T, S),
                     schedule = TimeInterval(save_interval),
                     filename = output_prefix * "_snapshots",
                     with_halos = true,
                     overwrite_existing = true)

# Let's goo!
@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

run!(simulation)

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
"""
