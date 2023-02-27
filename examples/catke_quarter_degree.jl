using ClimaOcean
using ClimaOcean.NearGlobalSimulations: quarter_degree_near_global_simulation, prettyelapsedtime
using ClimaOcean.DataWrangling: diffuse_tracers!

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: WallTimeInterval
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity

using JLD2
using CUDA
using Printf
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

#####
##### Boundary layer turbulence closure options
#####

# "Ri-based" --- uses calibrated defaults in Oceananigans
ri_based = RiBasedVerticalDiffusivity() 
default_catke = CATKEVerticalDiffusivity() 
clipped_catke = CATKEVerticalDiffusivity(maximum_diffusivity=1e-2) 
neutral_catke = ClimaOcean.neutral_catke
catke_and_convective_adj = (CATKEVerticalDiffusivity(), ConvectiveAdjustmentVerticalDiffusivity(convective_κz=1.0)) 
linear_equation_of_state = LinearEquationOfState(thermal_expansion=2e-4, haline_contraction=8e-5)
teos10 = TEOS10EquationOfState(; reference_density=1020)

# Choose closure
#vertical_mixing_closure = clipped_catke
#vertical_mixing_closure = default_catke
#vertical_mixing_closure = catke_and_convective_adj
vertical_mixing_closure = ri_based
equation_of_state = teos10 #linear_equation_of_state
arch = GPU()

dir = "/nobackup/users/glwagner/ClimaOcean"
#output_prefix = "default_catke"
output_prefix = "ri_based"
#output_prefix = "clipped_catke"
#output_prefix = "catke_flat_bottom"
#output_prefix = "catke_plus_conv_adj"

#####
##### Build the simulation
#####
                                                                                             
simulation = quarter_degree_near_global_simulation(arch;
                                                   stop_time = 10years,
                                                   minimum_ocean_depth = 10,
                                                   surface_temperature_relaxation_time_scale = 7days,
                                                   surface_salinity_relaxation_time_scale = 7days,
                                                   background_vertical_viscosity = 1e-4,
                                                   background_vertical_diffusivity = 1e-5,
                                                   equation_of_state,
                                                   vertical_mixing_closure)

T, S, e = simulation.model.tracers

@show T
@show S

Ty = compute!(Field(∂y(T)))
start = time_ns()
@info @sprintf("Starting max(∂y T) = %.2e", maximum(abs, Ty))

step(x, d, c) = 1/2 * (tanh((x - c) / d) + 1)
vertical_scale(x, y, z, t) = 10 + 190 * step(abs(z), 200, 1000)

diffuse_tracers!(simulation.model.grid;
                 tracers = (; T, S),
                 horizontal_scale = 100kilometers,
                 vertical_scale)

elapsed = 1e-9 * (time_ns() - start)
compute!(Ty)
@info @sprintf("Finished diffusing tracers (%s), max(∂y T) = %.2e", prettytime(elapsed), maximum(abs, Ty))

@show T
@show S

eᵢ(x, y, z) = 1e-6
set!(simulation.model, e=eᵢ)

# Define output
slices_save_interval = 1day
fields_save_interval = 30days
Nx, Ny, Nz = size(simulation.model.grid)

#=
simulation.output_writers[:checkpointer] = Checkpointer(simulation.model; dir,
                                                        prefix = output_prefix * "_checkpointer",
                                                        schedule = WallTimeInterval(60minutes),
                                                        cleanup = true,
                                                        overwrite_existing = true)
=#

model = simulation.model

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers); dir,
                                                      schedule = TimeInterval(slices_save_interval),
                                                      filename = output_prefix * "_fields",
                                                      with_halos = true,
                                                      overwrite_existing = true)

slice_indices = [(:, :, Nz), (:, :, Nz-20), (:, 300, :), (600, :, :)]
output_names = [:surface, :near_surface, :zonal, :meridional]
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
    outputs[:η] = Field(model.free_surface.η; indices)
    outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

    name = output_names[n]
    simulation.output_writers[name] = JLD2OutputWriter(model, outputs; dir,
                                                       #schedule = TimeInterval(slices_save_interval),
                                                       schedule = IterationInterval(10),
                                                       filename = output_prefix * "_fields_$name",
                                                       with_halos = true,
                                                       overwrite_existing = true)
end

@info "Running a $output_prefix simulation with Δt = $(prettytime(simulation.Δt))"

simulation.Δt = 1.0
simulation.stop_iteration = 1000
progress = simulation.callbacks[:progress].func
simulation.callbacks[:progress] = Callback(progress)
run!(simulation)

#=
simulation.Δt = 10.0
simulation.stop_iteration = 10000
run!(simulation)

simulation.Δt = 1minute
simulation.stop_iteration = 10000
run!(simulation)

simulation.Δt = 5minute
simulation.stop_iteration = Inf
run!(simulation)
=#

@info "Simulation took $(prettytime(simulation.run_wall_time))."

