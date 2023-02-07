using ClimaOcean.NearGlobalSimulations: quarter_degree_near_global_simulation

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

#####
##### Boundary layer turbulence closure options
#####

# "Ri-based" --- uses calibrated defaults in Oceananigans
ri_based = RiBasedVerticalDiffusivity() 
catke = CATKEVerticalDiffusivity() 

# Choose closure
boundary_layer_turbulence_closure = ri_based

#####
##### Build the simulation
#####

start_time = 0days
stop_time  = start_time + 10years

simulation = quarter_degree_near_global_simulation(GPU(); start_time, stop_time,
                                                   background_vertical_viscosity = 1e-2,
                                                   background_vertical_diffusivity = 1e-2,
                                                   boundary_layer_turbulence_closure = catke)

# Define output
slices_save_interval = 1day
fields_save_interval = 30days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "/nobackup/users/glwagner/ClimaOcean"
closure_name = typeof(boundary_layer_turbulence_closure).name.wrapper
output_prefix = "catke_near_global_$(Nx)_$(Ny)_$(Nz)"

simulation.output_writers[:checkpointer] = Checkpointer(simulation.model; dir,
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = WallTimeInterval(10minutes),
                                                        cleanup = true,
                                                        overwrite_existing = true)

model = simulation.model

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers); dir,
                                                      schedule = TimeInterval(slices_save_interval),
                                                      filename = output_prefix * "_fields",
                                                      with_halos = true,
                                                      overwrite_existing = true)

slice_indices = [(:, :, Nz), (:, :, Nz-20)]
output_names = [:surface, :near_surface]
for n = 1:2
    indices = slice_indices[n]

    outputs = Dict()

    for name in keys(model.tracers)
        c = model.tracers[name]
        outputs[name] = Field(c; indices)
    end

    outputs[:u] = Field(model.velocities.u; indices)
    outputs[:v] = Field(model.velocities.v; indices)
    outputs[:w] = Field(model.velocities.w; indices)
    outputs[:η] = model.free_surface.η
    outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

    name = output_names[n]
    simulation.output_writers[name] = JLD2OutputWriter(model, outputs; dir,
                                                       schedule = TimeInterval(slices_save_interval),
                                                       filename = output_prefix * "_fields_$name",
                                                       with_halos = true,
                                                       overwrite_existing = true)
end

@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"
simulation.Δt = 1minute
simulation.stop_iteration = 100
simulation.callbacks[:progress] = Callback(simulation.callbacks[:progress].func)
run!(simulation)

simulation.callbacks[:progress] = Callback(simulation.callbacks[:progress].func, IterationInterval(10))
simulation.stop_iteration = Inf
simulation.Δt = 10minutes
run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

