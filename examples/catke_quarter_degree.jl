using ClimaOcean
using ClimaOcean.NearGlobalSimulations: quarter_degree_near_global_simulation, prettyelapsedtime

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

using JLD2
using CUDA

#####
##### Boundary layer turbulence closure options
#####

# "Ri-based" --- uses calibrated defaults in Oceananigans
ri_based = RiBasedVerticalDiffusivity() 
default_catke = CATKEVerticalDiffusivity() 
neutral_catke = ClimaOcean.neutral_catke

# Choose closure
boundary_layer_turbulence_closure = ri_based

#####
##### Build the simulation
#####

simulation = quarter_degree_near_global_simulation(GPU();
                                                   stop_time = 10years,
                                                   background_horizontal_viscosity = 1e2,
                                                   background_horizontal_diffusivity = (T=0.0, S=0.0, e=1e2),
                                                   initial_vertical_diffusion_steps = 0,
                                                   initial_horizontal_diffusion_steps = 0,
                                                   background_vertical_viscosity = 1e-4,
                                                   background_vertical_diffusivity = 1e-5,
                                                   boundary_layer_turbulence_closure = default_catke)

#=
vitd = VerticallyImplicitTimeDiscretization()
vertical_diffusivity = VerticalScalarDiffusivity(vitd, ν=1e-4, κ=(T=1e-5, S=1e-5, e=1e-5))
horizontal_viscosity = HorizontalScalarDiffusivity(κ=(T=0.0, S=0.0, e=0.0), ν=1e3)
neutral_horizontal_diffusivity = HorizontalScalarDiffusivity(κ=(T=0.0, S=0.0, e=0.0), ν=0.0)

@info "Diffusing tracers..."; start=time_ns()
# Steps to diffuse tracer field
u, v, w = simulation.model.velocities
e = simulation.model.tracers.e
η = simulation.model.free_surface.η
simulation.Δt = 0.02
for step = 1:1000
    time_step!(simulation)
    fill!(parent(u), 0)
    fill!(parent(v), 0)
    fill!(parent(w), 0)
    fill!(parent(η), 0)
    fill!(parent(e), 0)
    @info "Step $step"
end
@info "    ... done diffusing tracers (" * prettyelapsedtime(start) * ")"

# Reset
simulation.model.clock.time = 0.0
simulation.model.clock.iteration = 0
simulation.model.closure = (default_catke, vertical_diffusivity, horizontal_viscosity)
=#

@info "Setting initial condition from file..."
filepath = "/nobackup/users/glwagner/ClimaOcean/near_global_1440_600_87_RiBasedVerticalDiffusivity_checkpointer.jld2"
file = jldopen(filepath)
u_data = file["u/data"]
v_data = file["v/data"]
T_data = file["T/data"]
S_data = file["S/data"]
η_data = file["η/data"]
time = file["clock"].time
close(file)

simulation.model.clock.time = time
parent(simulation.model.velocities.u) .= CuArray(u_data)
parent(simulation.model.velocities.v) .= CuArray(v_data)
parent(simulation.model.tracers.T) .= CuArray(T_data)
parent(simulation.model.tracers.S) .= CuArray(S_data)
parent(simulation.model.free_surface.η) .= CuArray(η_data)

eᵢ(x, y, z) = 1e-6
set!(simulation.model, e=eᵢ)

# Define output
slices_save_interval = 1day
fields_save_interval = 30days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "/nobackup/users/glwagner/ClimaOcean"
closure_name = typeof(boundary_layer_turbulence_closure).name.wrapper
output_prefix = "checkpoint_test_near_global_$(Nx)_$(Ny)_$(Nz)"

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

    u, v, w = model.velocities

    K = @at (Center, Center, Center) (u^2 + v^2) / 2
    outputs[:u] = Field(u; indices)
    outputs[:v] = Field(v; indices)
    outputs[:w] = Field(w; indices)
    outputs[:K] = Field(K; indices)
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
simulation.callbacks[:progress] = Callback(simulation.callbacks[:progress].func)
run!(simulation)

# simulation.callbacks[:progress] = Callback(simulation.callbacks[:progress].func, IterationInterval(10))
# simulation.stop_iteration = Inf
# simulation.Δt = 1minute
# run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

