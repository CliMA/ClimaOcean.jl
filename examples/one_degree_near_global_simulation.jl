using Oceananigans
using Oceananigans.Units
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation

# Build the simulation
simulation = one_degree_near_global_simulation()

# Define output
fields_save_interval = 5days
surface_fields_save_interval = 2hours
Nx, Ny, Nz = size(simulation.model.grid)
output_prefix = "near_global_$(Nx)_$(Ny)_$(Nz)"

model = simulation.model
outputs = merge(model.velocities, model.tracers)

simulation.output_writers[:fields] =
    JLD2OutputWriter(model, outputs,
                     schedule = TimeInterval(save_interval),
                     filename = output_prefix * "_fields",
                     with_halos = true,
                     overwrite_existing = true)

surface_outputs = Dict(name => Field(outputs[name], indices=(:, :, Nz)) for name in keys(outputs))
surface_outputs[:η] = model.free_surface.η
surface_outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices=(:, :, Nz))

simulation.output_writers[:surface_fields] =
    JLD2OutputWriter(model, surface_outputs,
                     schedule = TimeInterval(surface_fields_save_interval),
                     filename = output_prefix * "_surface_fields",
                     with_halos = true,
                     overwrite_existing = true)

@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"
run!(simulation)
@info "Simulation took $(prettytime(simulation.run_wall_time))."

