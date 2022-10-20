module NearGlobalSimulations

using Oceananigans
using Oceananigans.Units

using Oceananigans.Operators: Δzᵃᵃᶜ, ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.Coriolis: WetCellEnstrophyConservingScheme
using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity, FluxTapering

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity, MixingLength

using DataDeps
using Statistics
using JLD2
using Printf
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using CUDA: @allowscalar

import Oceananigans.Utils: prettysummary
import ..VerticalGrids

struct PiecewiseConstantVerticalDiffusivity <: Function
    z_transition :: Float64
    top_diffusivity :: Float64
    bottom_diffusivity :: Float64
end

@inline function (pcd::PiecewiseConstantVerticalDiffusivity)(x, y, z, t)
    zᵗ = pcd.z_transition
    κᵗ = pcd.top_diffusivity
    κᵇ = pcd.bottom_diffusivity
    return ifelse(z > zᵗ, κᵗ, κᵇ)
end

prettysummary(pcd::PiecewiseConstantVerticalDiffusivity) =
    string("PiecewiseConstantVerticalDiffusivity(",
           pcd.z_transition, ", ",
           pcd.top_diffusivity, ", ",
           pcd.bottom_diffusivity, ")")

Base.summary(pcd::PiecewiseConstantVerticalDiffusivity) = prettysummary(pcd)

const thirty_days = 30days

@inline current_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
@inline next_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
@inline cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / thirty_days, 1) * (u₂ - u₁)

@inline function surface_temperature_relaxation(i, j, grid, clock, fields, p)
    time = clock.time

    n₁ = current_time_index(time, p.Nmonths)
    n₂ = next_time_index(time, p.Nmonths)

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

    n₁ = current_time_index(time, p.Nmonths)
    n₂ = next_time_index(time, p.Nmonths)

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

@inline function surface_wind_stress(i, j, grid, clock, fields, p)
    time = clock.time
    n₁ = current_time_index(time, p.Nmonths)
    n₂ = next_time_index(time, p.Nmonths)

    @inbounds begin
        τ₁ = p.τ[i, j, n₁]
        τ₂ = p.τ[i, j, n₂]
    end

    return cyclic_interpolate(τ₁, τ₂, time)
end

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

"""
    one_degree_near_global_simulation(architecture = GPU(); kwargs...)

Return an Oceananigans.Simulation of Earth's ocean at 1 degree resolution.
"""
function one_degree_near_global_simulation(architecture = GPU();
    size                                         = (360, 150, 48),
    boundary_layer_turbulence_closure            = RiBasedVerticalDiffusivity(),
    background_vertical_diffusivity              = 1e-5,
    horizontal_viscosity                         = 5e4,
    surface_background_vertical_viscosity        = 1e-2,
    interior_background_vertical_viscosity       = 1e-4,
    vertical_viscosity_transition_depth          = 49.0,
    with_isopycnal_skew_symmetric_diffusivity    = true,
    surface_temperature_relaxation_time_scale    = 30days,
    surface_salinity_relaxation_time_scale       = 90days,
    isopycnal_κ_skew                             = 900.0,
    isopycnal_κ_symmetric                        = 900.0,
    max_isopycnal_slope                          = 1e-2,
    bottom_drag_coefficient                      = 3e-3,
    reference_density                            = 1029.0,
    reference_heat_capacity                      = 3991.0,
    reference_salinity                           = 34.0,
    time_step                                    = 20minutes,
    stop_iteration                               = Inf,
    start_time                                   = 345days,
    stop_time                                    = Inf,
    bathymetry_path                              = datadep"near_global_one_degree/bathymetry_lat_lon_360_150.jld2",
    initial_conditions_path                      = datadep"near_global_one_degree/initial_conditions_month_01_360_150_48.jld2",
    surface_boundary_conditions_path             = datadep"near_global_one_degree/surface_boundary_conditions_12_months_360_150.jld2",
    )

    size == (360, 150, 48) || throw(ArgumentError("Only size = (360, 150, 48) is supposed."))

    #####
    ##### Load surface boundary conditions and inital conditions
    ##### from ECCO version 4:
    ##### https://ecco.jpl.nasa.gov/drive/files
    #####
    ##### Bathymetry is interpolated from ETOPO1:
    ##### https://www.ngdc.noaa.gov/mgg/global/
    #####

    bathymetry_file = jldopen(bathymetry_path)
    bathymetry = bathymetry_file["bathymetry"]
    close(bathymetry_file)

    @info "Reading initial conditions..."; start=time_ns()
    initial_conditions_file = jldopen(initial_conditions_path)
    T_init = initial_conditions_file["T"]
    S_init = initial_conditions_file["S"]
    close(initial_conditions_file)
    @info "... read initial conditions (" * prettytime(1e-9 * (time_ns() - start)) * ")"

    # Files contain 12 arrays of monthly-averaged data from 1992
    @info "Reading boundary conditions..."; start=time_ns()
    boundary_conditions_file = jldopen(surface_boundary_conditions_path)
    τˣ = - boundary_conditions_file["τˣ"] ./ reference_density
    τʸ = - boundary_conditions_file["τʸ"] ./ reference_density
    T★ = + boundary_conditions_file["Tₛ"]
    S★ = + boundary_conditions_file["Sₛ"]
    Q★ = - boundary_conditions_file["Qᶠ"] ./ reference_density ./ reference_heat_capacity
    F★ = - boundary_conditions_file["Sᶠ"] ./ reference_density .* reference_salinity
    close(boundary_conditions_file)
    @info "... read boundary conditions (" * prettytime(1e-9 * (time_ns() - start)) * ")"

    # Convert boundary conditions arrays to GPU
    τˣ = arch_array(architecture, τˣ)
    τʸ = arch_array(architecture, τʸ)
    target_sea_surface_temperature = T★ = arch_array(architecture, T★)
    target_sea_surface_salinity    = S★ = arch_array(architecture, S★)
    surface_temperature_flux       = Q★ = arch_array(architecture, Q★)
    surface_salt_flux              = F★ = arch_array(architecture, F★)

    # Stretched faces from ECCO Version 4 (49 levels in the vertical)
    z_faces = VerticalGrids.z_49_levels_10_to_400_meter_spacing

    # A spherical domain
    underlying_grid = LatitudeLongitudeGrid(architecture; size,
                                            longitude = (-180, 180),
                                            latitude = (-75, 75),
                                            halo = (5, 5, 5),
                                            z = z_faces)

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))

    @info "Created $grid"

    #####
    ##### Physics and model setup
    #####

    νz = PiecewiseConstantVerticalDiffusivity(-vertical_viscosity_transition_depth,
                                              surface_background_vertical_viscosity,
                                              interior_background_vertical_viscosity)

    vitd = VerticallyImplicitTimeDiscretization()

    horizontal_viscosity = HorizontalScalarDiffusivity(ν=horizontal_viscosity)
    vertical_viscosity   = VerticalScalarDiffusivity(vitd, ν=νz, κ=background_vertical_diffusivity)

    closures = Any[horizontal_viscosity, boundary_layer_turbulence_closure, vertical_viscosity]

    if with_isopycnal_skew_symmetric_diffusivity
        issd = IsopycnalSkewSymmetricDiffusivity(κ_skew = isopycnal_κ_skew,
                                                 κ_symmetric = isopycnal_κ_symmetric,
                                                 slope_limiter = FluxTapering(max_isopycnal_slope))
        push!(closures, issd)
    end

    # TODO: do this internally in model constructor
    closures = tuple(closures...)

    #####
    ##### Boundary conditions / time-dependent fluxes 
    #####

    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = bottom_drag_coefficient)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = bottom_drag_coefficient)

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

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

    Nmonths = 12 # number of months in the forcing file
    u_wind_stress_parameters = (; τ=τˣ, Nmonths)
    v_wind_stress_parameters = (; τ=τʸ, Nmonths)
    u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form=true, parameters=u_wind_stress_parameters)
    v_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form=true, parameters=v_wind_stress_parameters)

    Δz_top = @allowscalar Δzᵃᵃᶜ(1, 1, grid.Nz, grid.underlying_grid)

    T_relaxation_parameters = (; λ = Δz_top / surface_temperature_relaxation_time_scale,
                                 Nmonths,
                                 T★ = target_sea_surface_temperature,
                                 Q★ = surface_temperature_flux)

    S_relaxation_parameters = (; λ = Δz_top / surface_salinity_relaxation_time_scale,
                                 Nmonths,
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

    # Because MITgcm forcing starts at Jan 15 (?)
    model.clock.time = start_time

    simulation = Simulation(model; Δt=time_step, stop_iteration, stop_time)

    start_time = [time_ns()]

    function progress(sim)
        wall_time = (time_ns() - start_time[1]) * 1e-9

        u = sim.model.velocities.u
        w = sim.model.velocities.w

        intw  = Array(interior(w))
        max_w = findmax(intw)

        mw = max_w[1]
        iw = max_w[2]

        msg1 = @sprintf("Time: % 12s, iteration: %d, ", prettytime(sim), iteration(sim))
        msg2 = @sprintf("max(|u|): %.2e ms⁻¹, wmax: %.2e, loc: (%d, %d, %d), ",
                        maximum(abs, u), mw, iw[1], iw[2], iw[3])
        msg3 = @sprintf("wall time: %s", prettytime(wall_time))

        @info msg1 * msg2 * msg3

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    return simulation
end

end # module

