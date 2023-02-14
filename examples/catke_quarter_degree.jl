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
using Printf

#####
##### Boundary layer turbulence closure options
#####

# "Ri-based" --- uses calibrated defaults in Oceananigans
ri_based = RiBasedVerticalDiffusivity() 
default_catke = CATKEVerticalDiffusivity() 
neutral_catke = ClimaOcean.neutral_catke

# Choose closure
boundary_layer_turbulence_closure = ri_based

arch = GPU()

#####
##### Build the simulation
#####

linear_equation_of_state = LinearEquationOfState(thermal_expansion=2e-4, haline_contraction=8e-5)
                                                                                             
simulation = quarter_degree_near_global_simulation(arch;
                                                   stop_time = 10years,
                                                   background_vertical_viscosity = 1e-4,
                                                   background_vertical_diffusivity = 1e-5,
                                                   boundary_layer_turbulence_closure = ri_based)

function diffuse_tracers!(grid;
                          tracers,
                          horizontal_scale = 0,
                          vertical_scale = 0,
                          fractional_time_step = 1e-2)

    ϵ = fractional_time_step
    κh = horizontal_scale^2
    κz = vertical_scale^2

    grid = simulation.model.grid
    Nx, Ny, Nz = size(grid)
    Az = minimum(grid.Azᶜᶜᵃ[1:Ny])
    Δt = ϵ * Az / κh
    @show Nt = ceil(Int, 1 / Δt)

    vitd = VerticallyImplicitTimeDiscretization()
    vertical_smoothing = VerticalScalarDiffusivity(vitd, κ = κz^2)
    horizontal_smoothing = HorizontalScalarDiffusivity(κ = κh^2)

    smoothing_model = HydrostaticFreeSurfaceModel(; grid, tracers,
                                                  velocities = PrescribedVelocityFields(),
                                                  buoyancy = nothing,
                                                  momentum_advection = nothing,
                                                  tracer_advection = nothing,
                                                  coriolis = nothing,
                                                  closure = (horizontal_smoothing, vertical_smoothing))

    @info string("Smoothing tracers ", keys(tracers))

    for n = 1:Nt
        time_step!(smoothing_model, Δt)
    end

    return nothing
end

grid = simulation.model.grid
T, S, e = simulation.model.tracers
T_dummy = CenterField(grid, data=T.data)
S_dummy = CenterField(grid, data=S.data)

Ty = compute!(Field(∂y(T)))
@info @sprintf("Starting max(∂y T) = %.2e", maximum(abs, Ty))

start = time_ns()

diffuse_tracers!(simulation.model.grid;
                 tracers = (T=T_dummy, S=S_dummy),
                 horizontal_scale = 100kilometers,
                 vertical_scale = 10meters)

elapsed = 1e-9 * (time_ns() - start)
compute!(Ty)

@info @sprintf("Finished diffusing tracers (%s), max(∂y T) = %.2e",
               prettytime(elapsed), maximum(abs, Ty))

eᵢ(x, y, z) = 1e-6
set!(simulation.model, e=eᵢ)

# Define output
slices_save_interval = 1day
fields_save_interval = 30days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "/nobackup/users/glwagner/ClimaOcean"
closure_name = typeof(boundary_layer_turbulence_closure).name.wrapper
#output_prefix = "catke_test_near_global_$(Nx)_$(Ny)_$(Nz)"
output_prefix = "ri_based_test_near_global_$(Nx)_$(Ny)_$(Nz)"

simulation.output_writers[:checkpointer] = Checkpointer(simulation.model; dir,
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = WallTimeInterval(20minutes),
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
                                                       #schedule = TimeInterval(slices_save_interval),
                                                       schedule = IterationInterval(1),
                                                       filename = output_prefix * "_fields_$name",
                                                       with_halos = true,
                                                       overwrite_existing = true)
end

@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

simulation.Δt = 10minutes
simulation.stop_iteration = 100
simulation.callbacks[:progress] = Callback(simulation.callbacks[:progress].func)
run!(simulation)

# simulation.callbacks[:progress] = Callback(simulation.callbacks[:progress].func, IterationInterval(10))
# simulation.stop_iteration = Inf
# simulation.Δt = 1minute
# run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."
