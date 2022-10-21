using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation

# Build the simulation
with_isopycnal_skew_symmetric_diffusivity = true
start_time = 345days
simulation = one_degree_near_global_simulation(; start_time, with_isopycnal_skew_symmetric_diffusivity, 
                                               stop_time = start_time+10years)

# Define output
slices_save_interval = 5days
fields_save_interval = 10days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "./" 
output_prefix = "near_global_$(Nx)_$(Ny)_$(Nz)"
with_isopycnal_skew_symmetric_diffusivity || (output_prefix *= "_no_gm")

simulation.output_writers[:checkpointer] = Checkpointer(simulation.model; dir,
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = WallTimeInterval(10minutes),
                                                        cleanup = true,
                                                        overwrite_existing = true)

model = simulation.model

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers); dir,
                                                      schedule = TimeInterval(fields_save_interval),
                                                      filename = output_prefix * "_fields",
                                                      with_halos = true,
                                                      overwrite_existing = true)

surface_outputs = Dict()
indices = (:, :, Nz)
surface_outputs[:u] = Field(model.velocities.u; indices)
surface_outputs[:v] = Field(model.velocities.v; indices)
surface_outputs[:w] = Field(model.velocities.w; indices)
surface_outputs[:T] = Field(model.tracers.T; indices)
surface_outputs[:S] = Field(model.tracers.S; indices)
surface_outputs[:η] = model.free_surface.η
surface_outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

simulation.output_writers[:surface] = JLD2OutputWriter(model, surface_outputs; dir,
                                                       schedule = TimeInterval(slices_save_interval),
                                                       filename = output_prefix * "_surface_fields",
                                                       with_halos = true,
                                                       overwrite_existing = true)

#=
# Add vertical slices or "transects"
transects = (pacific = (11, :, :),
             atlantic = (153, :, :),
             circumantarctic = (:, 16, :))

b = buoyancy(model)
outputs = merge(model.velocities, model.tracers, (; b))

for name in keys(transects)
    indices = transects[name]
    outputs = Dict{Symbol, Any}(n => Field(outputs[n]; indices) for n in keys(outputs))
    outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

    simulation.output_writers[:name] = JLD2OutputWriter(model, outputs; dir,
                                                        schedule = TimeInterval(slices_save_interval),
                                                        filename = output_prefix * "_$name",
                                                        with_halos = true,
                                                        overwrite_existing = true)
end
=#

@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

