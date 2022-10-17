using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation

# Build the simulation
simulation = one_degree_near_global_simulation()
simulation.stop_iteration = Inf
simulation.stop_time = simulation.model.clock.time + 10years

# Define output
save_interval = 5days
Nx, Ny, Nz = size(simulation.model.grid)
output_prefix = "near_global_$(Nx)_$(Ny)_$(Nz)"

simulation.output_writers[:checkpointer] = Checkpointer(simulation.model,
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = WallTimeInterval(10minutes),
                                                        cleanup = true,
                                                        overwrite_existing = true)

model = simulation.model
surface_outputs = Dict()
indices = (:, :, Nz)
surface_outputs[:u] = Field(model.velocities.u; indices)
surface_outputs[:v] = Field(model.velocities.v; indices)
surface_outputs[:w] = Field(model.velocities.w; indices)
surface_outputs[:T] = Field(model.tracers.T; indices)
surface_outputs[:S] = Field(model.tracers.S; indices)
surface_outputs[:η] = model.free_surface.η
surface_outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

simulation.output_writers[:surface] = JLD2OutputWriter(model, surface_outputs,
                                                       schedule = TimeInterval(save_interval),
                                                       filename = output_prefix * "_surface_fields",
                                                       with_halos = true,
                                                       overwrite_existing = true)

# Add vertical slices or "transects"
transects = (pacific = (11, :, :),
             atlantic = (153, :, :),
             circumantarctic = (:, 16, :))

outputs = merge(model.velocities, model.tracers)

for name in keys(transects)
    indices = transects[name]
    outputs = Dict(n => Field(outputs[n]; indices) for n in keys(outputs))

    simulation.output_writers[:name] = JLD2OutputWriter(model, outputs,
                                                        schedule = TimeInterval(save_interval),
                                                        filename = output_prefix * "_$name",
                                                        with_halos = true,
                                                        overwrite_existing = true)
end

@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

