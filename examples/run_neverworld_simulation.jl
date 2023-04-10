using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ClimaOcean.IdealizedSimulations: neverworld_simulation
using Printf

minimum_turbulent_kinetic_energy = 1e-6
minimum_convective_buoyancy_flux = 1e-11
closure = CATKEVerticalDiffusivity(; minimum_turbulent_kinetic_energy,
                                   minimum_convective_buoyancy_flux)

simulation = neverworld_simulation(GPU();
                                   horizontal_size = (240, 280),
                                   longitude = (0, 60),
                                   latitude = (-70, 0),
                                   time_step = 10minutes,
                                   stop_time = 200 * 360days,
                                   closure)

model = simulation.model
grid = model.grid

@show grid

start_time = Ref(time_ns())

function progress(sim) 
    b = sim.model.tracers.b
    e = sim.model.tracers.e
    u, v, w = sim.model.velocities

    msg = @sprintf("Iter: %d, time: %s, extrema(b): (%6.2e, %6.2e)",
                   iteration(sim), prettytime(sim), minimum(b), maximum(b))

    msg *= @sprintf(", extrema(e): (%6.2e, %6.2e)", minimum(e), maximum(e))

    msg *= @sprintf(", max|u|: %6.2e, max|w|: %6.2e",
                    maximum(maximum(abs, q) for q in (u, v, w)), maximum(abs, w))

    try 
        κᶜ = sim.model.diffusivity_fields.Kᶜ
        msg *= @sprintf(", extrema(κᶜ): (%6.2e, %6.2e)", minimum(κᶜ), maximum(κᶜ))
    catch
    end

    elapsed = 1e-9 * (time_ns() - start_time[])
    msg *= @sprintf(", wall time: %s", prettytime(elapsed))
    start_time[] = time_ns()

    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Set up output
Nx, Ny, Nz = size(grid)
i = round(Int, Nx/10) # index for yz-sliced output

diffusivity_fields = (; κᶜ = model.diffusivity_fields.κᶜ)
outputs = merge(model.velocities, model.tracers, diffusivity_fields)
zonally_averaged_outputs = NamedTuple(n => Average(outputs[n], dims=1) for n in keys(outputs))

simulation.output_writers[:yz] = JLD2OutputWriter(model, outputs;
                                                  schedule = TimeInterval(6days),
                                                  filename = "neverworld_yz.jld2",
                                                  indices = (i, :, :),
                                                  with_halos = true,
                                                  overwrite_existing = true)

simulation.output_writers[:zonal] = JLD2OutputWriter(model, zonally_averaged_outputs;
                                                     schedule = TimeInterval(6days),
                                                     filename = "neverworld_zonal_average.jld2",
                                                     with_halos = true,
                                                     overwrite_existing = true)

simulation.output_writers[:xy] = JLD2OutputWriter(model, outputs;
                                                  schedule = TimeInterval(6days),
                                                  filename = "neverworld_xy.jld2",
                                                  indices = (:, :, Nz),
                                                  with_halos = true,
                                                  overwrite_existing = true)

simulation.output_writers[:xyz] = JLD2OutputWriter(model, outputs;
                                                   schedule = TimeInterval(90days),
                                                   filename = "neverworld_xyz.jld2",
                                                   with_halos = true,
                                                   overwrite_existing = true)

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(360days),
                                                        prefix = "neverworld_checkpoint",
                                                        cleanup = true)

run!(simulation)

