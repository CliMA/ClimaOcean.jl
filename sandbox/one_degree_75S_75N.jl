using ClimaOcean

using Oceananigans
using Oceananigans.Units

using Oceananigans.Operators: Δzᵃᵃᶜ
using Oceananigans.Architectures: on_architecture
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.Coriolis: WetCellEnstrophyConservingScheme
using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity, FluxTapering

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity, MixingLength

using Statistics
using JLD2
using Printf
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using CUDA: @allowscalar

# 1 degree resolution
Nx = 360
Ny = 150
Nz = 48

arch = GPU()
reference_density        = 1029.0
reference_heat_capacity  = 3991.0
reference_salinity       = 34.0

output_prefix = "./near_global_lat_lon_$(Nx)_$(Ny)_$(Nz)_fine"
pickup = false
save_interval = 5days

#####
##### Load surface boundary conditions and inital conditions
##### from ECCO version 4:
##### https://ecco.jpl.nasa.gov/drive/files
#####
##### Bathymetry is interpolated from ETOPO1:
##### https://www.ngdc.noaa.gov/mgg/global/
#####

bathymetry_file = jldopen("bathymetry_360_150_75S_75N.jld2")
bathymetry = bathymetry_file["bathymetry"]
close(bathymetry_file)

@info "Reading initial conditions..."; start=time_ns()
initial_conditions_file = jldopen("initial_conditions_360_150_48_75S_75N.jld2")
T_init = initial_conditions_file["T"]
S_init = initial_conditions_file["S"]
close(initial_conditions_file)
@info "... read initial conditions (" * prettytime(1e-9 * (time_ns() - start)) * ")"

# Files contain 12 arrays of monthly-averaged data from 1992
@info "Reading boundary conditions..."; start=time_ns()
boundary_conditions_file = jldopen("surface_boundary_conditions_360_150_75S_75N.jld2")
τˣ = - boundary_conditions_file["τˣ"] ./ reference_density
τʸ = - boundary_conditions_file["τʸ"] ./ reference_density
T★ = + boundary_conditions_file["Tₛ"]
S★ = + boundary_conditions_file["Sₛ"]
Q★ = - boundary_conditions_file["Qᶠ"] ./ reference_density ./ reference_heat_capacity
F★ = - boundary_conditions_file["Sᶠ"] ./ reference_density .* reference_salinity
close(boundary_conditions_file)
@info "... read boundary conditions (" * prettytime(1e-9 * (time_ns() - start)) * ")"

# Convert boundary conditions arrays to GPU
τˣ = on_architecture(arch, τˣ)
τʸ = on_architecture(arch, τʸ)
target_sea_surface_temperature = T★ = on_architecture(arch, T★)
target_sea_surface_salinity = S★ = on_architecture(arch, S★)
surface_temperature_flux = Q★ = on_architecture(arch, Q★)
surface_salt_flux = F★ = on_architecture(arch, F★)

# Stretched faces from ECCO Version 4 (49 levels in the vertical)
z_faces = ClimaOcean.VerticalGrids.z_49_levels_10_to_400_meter_spacing

# A spherical domain
underlying_grid = LatitudeLongitudeGrid(arch,
                                        size = (Nx, Ny, Nz),
                                        longitude = (-180, 180),
                                        latitude = (-75, 75),
                                        halo = (5, 5, 5),
                                        z = z_faces)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))

@info "Created $grid"

#####
##### Physics and model setup
#####

const surface_νz = 1e-2
const background_νz = 1e-4
const background_κz = 1e-5

@inline νz(x, y, z, t) = ifelse(z > -49, surface_νz, background_νz)

horizontal_viscosity = HorizontalScalarDiffusivity(ν=5e4)
vertical_mixing      = RiBasedVerticalDiffusivity()
vertical_viscosity   = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(),
                                                 ν=νz, κ=background_κz)

κ_skew = 900.0      # [m² s⁻¹] skew diffusivity
κ_symmetric = 900.0 # [m² s⁻¹] symmetric diffusivity

gent_mcwilliams_diffusivity = IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric,
                                                                slope_limiter = FluxTapering(1e-2))

#closures = gent_mcwilliams_diffusivity
    
closures = (
    vertical_viscosity,
    horizontal_viscosity,
    vertical_mixing,
    gent_mcwilliams_diffusivity,
)

#####
##### Boundary conditions / time-dependent fluxes 
#####

const Nyears = 2.0
const Nmonths = 12
const thirty_days = 30days

@inline current_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
@inline next_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
@inline cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / thirty_days, 1) * (u₂ - u₁)

Δz_top = @allowscalar Δzᵃᵃᶜ(1, 1, grid.Nz, grid.underlying_grid)
Δz_bottom = @allowscalar Δzᵃᵃᶜ(1, 1, 1, grid.underlying_grid)

@inline function surface_wind_stress(i, j, grid, clock, fields, τ)
    time = clock.time
    n₁ = current_time_index(time, Nmonths)
    n₂ = next_time_index(time, Nmonths)

    @inbounds begin
        τ₁ = τ[i, j, n₁]
        τ₂ = τ[i, j, n₂]
    end

    return cyclic_interpolate(τ₁, τ₂, time)
end

Δz_top = @allowscalar grid.Δzᵃᵃᶜ[Nz]

using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ

# Linear bottom drag:
μ = 0.003 # Non dimensional

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = @inbounds sqrt(fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)^2)
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = @inbounds sqrt(fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)^2)

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields)

drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ)
drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ)

no_slip_bc = ValueBoundaryCondition(0)

u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u,
                                          west = no_slip_bc,
                                          east = no_slip_bc,
                                          south = no_slip_bc,
                                          north = no_slip_bc)

v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v,
                                          west = no_slip_bc,
                                          east = no_slip_bc,
                                          south = no_slip_bc,
                                          north = no_slip_bc)

u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ)
v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ)

u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τˣ);
v_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τʸ);

@inline function surface_temperature_relaxation(i, j, grid, clock, fields, p)
    time = clock.time

    n₁ = current_time_index(time, Nmonths)
    n₂ = next_time_index(time, Nmonths)

    @inbounds begin
        T★₁ = p.T★[i, j, n₁]
        T★₂ = p.T★[i, j, n₂]
        Q★₁ = p.Q★[i, j, n₁]
        Q★₂ = p.Q★[i, j, n₂]
        T_surface = fields.T[i, j, grid.Nz]
    end

    T★ = cyclic_interpolate(T★₁, T★₂, time)
    Q★ = cyclic_interpolate(Q★₁, Q★₂, time)

    return Q★ + p.λ * (T_surface - T★)
end

@inline function surface_salinity_relaxation(i, j, grid, clock, fields, p)
    time = clock.time

    n₁ = current_time_index(time, Nmonths)
    n₂ = next_time_index(time, Nmonths)

    @inbounds begin
        S★₁ = p.S★[i, j, n₁]
        S★₂ = p.S★[i, j, n₂]
        F★₁ = p.F★[i, j, n₁]
        F★₂ = p.F★[i, j, n₂]
        S_surface = fields.S[i, j, grid.Nz]
    end

    S★ = cyclic_interpolate(S★₁, S★₂, time)
    F★ = cyclic_interpolate(F★₁, F★₂, time)

    return - F★ + p.λ * (S_surface - S★)
end

T_relaxation_parameters = (λ = Δz_top/30days,
                           T★ = target_sea_surface_temperature,
                           Q★ = surface_temperature_flux)

S_relaxation_parameters = (λ = Δz_top/90days,
                           S★ = target_sea_surface_salinity,
                           F★ = surface_salt_flux)

T_surface_relaxation_bc = FluxBoundaryCondition(surface_temperature_relaxation,
                                                discrete_form = true,
                                                parameters = T_relaxation_parameters)

S_surface_relaxation_bc = FluxBoundaryCondition(surface_salinity_relaxation,
                                                discrete_form = true,
                                                parameters = S_relaxation_parameters)

u_bcs = FieldBoundaryConditions(top = u_wind_stress_bc,
                                bottom = u_bottom_drag_bc,
                                immersed = u_immersed_bc)

v_bcs = FieldBoundaryConditions(top = v_wind_stress_bc,
                                bottom = v_bottom_drag_bc,
                                immersed = v_immersed_bc)

T_bcs = FieldBoundaryConditions(top = T_surface_relaxation_bc)
S_bcs = FieldBoundaryConditions(top = S_surface_relaxation_bc)

equation_of_state = TEOS10EquationOfState(; reference_density)
#equation_of_state = LinearEquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)
coriolis = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme())
free_surface = ImplicitFreeSurface()

@info "Building a model..."; start=time_ns()

model = HydrostaticFreeSurfaceModel(; grid, free_surface, buoyancy, coriolis,
                                    momentum_advection = VectorInvariant(), 
                                    tracer_advection = WENO(underlying_grid),
                                    closure = closures,
                                    boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
                                    tracers = (:T, :S))
@info "... built $model."
@info "Model building time: " * prettytime(1e-9 * (time_ns() - start))

#####
##### Initial condition:
#####

set!(model, T=T_init, S=S_init)

# Because JMC's forcing starts at Jan 15
model.clock.time = 345days

#####
##### Simulation setup
#####

Δt = 20minutes
simulation = Simulation(model; Δt, stop_iteration=100) #stop_time=Nyears * years)

start_time = [time_ns()]

function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    u = sim.model.velocities.u
    w = sim.model.velocities.w

    intw  = Array(interior(w))
    max_w = findmax(intw)

    mw = max_w[1]
    iw = max_w[2]

    @info @sprintf("Time: % 12s, iteration: %d, max(|u|): %.2e ms⁻¹, wmax: %.2e , loc: (%d, %d, %d), wall time: %s", 
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration, maximum(abs, u), mw, iw[1], iw[2], iw[3], 
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S

simulation.output_writers[:surface_fields] =
    JLD2OutputWriter(model, (; u, v, w, T, S),
                     schedule = TimeInterval(save_interval),
                     filename = output_prefix * "_snapshots",
                     with_halos = true,
                     overwrite_existing = true)

# Let's goo!
@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

run!(simulation; pickup)

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""
