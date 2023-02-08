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
                                                   initial_vertical_diffusion_steps = 100,
                                                   initial_horizontal_diffusion_steps = 10000,
                                                   boundary_layer_turbulence_closure = ri_based)

# Define output
slices_save_interval = 1day
fields_save_interval = 30days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "/nobackup/users/glwagner/ClimaOcean"
closure_name = typeof(boundary_layer_turbulence_closure).name.wrapper
output_prefix = "smoothed_ic_near_global_$(Nx)_$(Ny)_$(Nz)_$closure_name"
pickup = false
overwrite_existing = !pickup

simulation.output_writers[:checkpointer] = Checkpointer(simulation.model; dir,
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = WallTimeInterval(10minutes),
                                                        cleanup = true,
                                                        overwrite_existing)

model = simulation.model

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers); dir,
                                                      schedule = TimeInterval(fields_save_interval),
                                                      filename = output_prefix * "_fields",
                                                      with_halos = true,
                                                      overwrite_existing)

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
                                                       overwrite_existing)
end

@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"
simulation.Δt = 5minute
run!(simulation; pickup)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

