using ClimaOcean.NearGlobalSimulations: one_degree_near_global_simulation

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    MixingLength, TurbulentKineticEnergyEquation, CATKEVerticalDiffusivity

using ParameterEstimocean.Parameters: closure_with_parameters
using JLD2

#####
##### Boundary layer turbulence closure options
#####

# "Ri-based" --- uses calibrated defaults in Oceananigans
ri_based = RiBasedVerticalDiffusivity() 
default_catke = CATKEVerticalDiffusivity()

# Choose closure
boundary_layer_turbulence_closure = ri_based

@show boundary_layer_turbulence_closure

#####
##### Build the simulation
#####

start_time = 0days
stop_time = start_time + 10years
with_isopycnal_skew_symmetric_diffusivity = true

simulation = one_degree_near_global_simulation(; start_time, stop_time,
    with_isopycnal_skew_symmetric_diffusivity,
    boundary_layer_turbulence_closure,
    isopycnal_κ_skew = 900.0,
    isopycnal_κ_symmetric = 900.0,
    interior_background_vertical_viscosity = 1e-4,
    surface_background_vertical_viscosity = 1e-4,
)

# Define output
slices_save_interval = 2day
fields_save_interval = 10days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "/nobackup/users/glwagner/ClimaOcean"
closure_name = typeof(boundary_layer_turbulence_closure).name.wrapper
output_prefix = "near_global_$(Nx)_$(Ny)_$(Nz)_$closure_name"
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

slice_indices = [(:, :, Nz), (:, :, Nz-10)]
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

simulation.Δt = 1minute
simulation.stop_time = time(simulation) + 2days
run!(simulation)

simulation.stop_time = start_time + 1.2years
simulation.Δt = 20minutes
run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

