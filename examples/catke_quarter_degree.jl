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

arch = GPU()

#####
##### Build the simulation
#####

simulation = quarter_degree_near_global_simulation(arch;
                                                   stop_time = 10years,
                                                   background_vertical_viscosity = 1e-4,
                                                   background_vertical_diffusivity = 1e-5,
                                                   boundary_layer_turbulence_closure = default_catke)

grid = simulation.model.grid
Nx, Ny, Nz = size(grid)
Δh = 200kilometers
Az = minimum(grid.Azᶜᶜᵃ[1:Ny])
Δz = 10meters
ϵt = 5e-2
@show Δt = ϵt * Az / Δh^2
Nt = ceil(Int, Δh / sqrt(Az) / ϵt)

vitd = VerticallyImplicitTimeDiscretization()
vertical_smoothing = VerticalScalarDiffusivity(vitd, κ = Δz^2)
horizontal_smoothing = HorizontalScalarDiffusivity(κ = Δh^2)

smoothing_model = HydrostaticFreeSurfaceModel(grid = simulation.model.grid,
                                              #velocities = PrescribedVelocityFields(),
                                              buoyancy = nothing,
                                              momentum_advection = nothing,
                                              tracer_advection = nothing,
                                              coriolis = nothing,
                                              tracers = simulation.model.tracers,
                                              closure = (horizontal_smoothing, vertical_smoothing))

T = simulation.model.tracers.T
S = simulation.model.tracers.S

@info "Smoothing initial condition..."
for n = 1:Nt
    start = time_ns()
    time_step!(smoothing_model, Δt)
    elapsed = 1e-9 * (time_ns() - start)
    @info "Step $n of $Nt, " * prettytime(elapsed)

    @info maximum(abs, ∂y(T))
end

Ty = compute!(Field(∂y(T)))
Sy = compute!(Field(∂y(S)))

Ty = compute!(Field(∂y(T)))
Sy = compute!(Field(∂y(S)))

Nx, Ny, Nz = size(T)

using Oceananigans.ImmersedBoundaries: mask_immersed_field!

mask_immersed_field!(T, NaN)
mask_immersed_field!(S, NaN)
mask_immersed_field!(Ty, NaN)
mask_immersed_field!(Sy, NaN)

#=
using GLMakie

fig = Figure(resolution=(1800, 1200))
axT = Axis(fig[1, 1])
axS = Axis(fig[2, 1])

axTy = Axis(fig[1, 2])
axSy = Axis(fig[2, 2])

heatmap!(axT, interior(T, :, :, Nz))
heatmap!(axS, interior(S, :, :, Nz))

mask_immersed_field!(Ty, 0)
mask_immersed_field!(Sy, 0)
Tylim = maximum(abs, Ty) / 100
Sylim = maximum(abs, Sy) / 100

heatmap!(axTy, interior(Ty, :, :, Nz), colorrange=(-Tylim, Tylim), colormap=:balance)
heatmap!(axSy, interior(Sy, :, :, Nz), colorrange=(-Sylim, Sylim), colormap=:balance)

display(fig)
=#

eᵢ(x, y, z) = 1e-6
set!(simulation.model, e=eᵢ)

# Define output
slices_save_interval = 1day
fields_save_interval = 30days
Nx, Ny, Nz = size(simulation.model.grid)

dir = "/nobackup/users/glwagner/ClimaOcean"
closure_name = typeof(boundary_layer_turbulence_closure).name.wrapper
output_prefix = "catke_test_near_global_$(Nx)_$(Ny)_$(Nz)"

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
=#
