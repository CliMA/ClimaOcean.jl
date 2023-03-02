using Oceananigans.Grids: ynode
using Oceananigans.Operators: Δzᵃᵃᶜ

# Slightly off-center vortical perturbations
ψ̃(x, y, ℓ, k) = exp(-(y + ℓ/10)^2 / 2ℓ^2) * cos(k * x) * cos(k * y)

# Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
ũ(x, y, ℓ, k) = + ψ̃(x, y, ℓ, k) * (k * tan(k * y) + y / ℓ^2)
ṽ(x, y, ℓ, k) = - ψ̃(x, y, ℓ, k) * k * tan(k * x)

@inline function aqua_planet_wind_stress(i, j, grid, clock, fields, parameters)
    φ = ynode(Center(), j, grid)
    τ₀ = parameters.maximum_wind_stress
    Δφₑ = parameters.equatorial_wind_extent
    Δφ = parameters.planet_extent
    return τ₀ * (cos(3π * φ/Δφ) - exp(-φ^2 / 2Δφₑ^2))
end

@inline function aqua_planet_sea_surface_temperature(φ, parameters)
    Tₑ = parameters.equilibrium_equatorial_temperature
    Tₚ = parameters.equilibrium_polar_temperature
    Δφ = parameters.planet_extent
    ϵ = parameters.temperature_asymmetry_factor
    φ₀ = parameters.temperature_asymmetry_latitude
    return Tₚ + (Tₑ - Tₚ) * (cos(π * φ/Δφ)^2 + ϵ * exp(-(φ - φ₀)^2 / 2φ^2))
end

struct AquaPlanetInitialCondition{P} <: Function
    scale_height :: Float64
    bottom_temperature :: Float64
    sea_surface_temperature_parameters :: P
end

function (ic::AquaPlanetInitialCondition)(λ, φ, z)
    parameters = ic.sea_surface_temperature_parameters
    T⁰ = aqua_planet_sea_surface_temperature(φ, parameters)
    Tᵇ = min(ic.bottom_temperature, T⁰) # prevent unstable stratification
    return Tᵇ + (T⁰ - Tᵇ) * exp(z / ic.scale_height)
end

@inline function aqua_planet_temperature_flux(i, j, grid, clock, fields, parameters)
    k = grid.Nz
    T⁰ = @inbounds fields.T[i, j, k]
    φ = ynode(Center(), j, grid)
    T★ = aqua_planet_sea_surface_temperature(φ, parameters)
    λ = parameters.surface_temperature_relaxation_time_scale
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
    q★ = Δz / λ
    return q★ * (T⁰ - T★)
end

"""
    near_global_aqua_planet_simulation(architecture = GPU(); kwargs...)

Return an Oceananigans.Simulation of an aqua planet.
"""
function aqua_planet_simulation(architecture = GPU();
    size                                         = (90, 30, 24),
    longitude                                    = (-180, 180),
    latitude                                     = (-60, 60),
    z                                            = (-3000, 0),
    closure                                      = RiBasedVerticalDiffusivity(),
    surface_temperature_relaxation_time_scale    = 10days,
    surface_salinity_relaxation_time_scale       = 90days,
    maximum_wind_stress                          = 1e-4,
    equatorial_wind_extent                       = 10.0,  # degrees latitude
    temperature_asymmetry_factor                 = 0.1,
    temperature_asymmetry_latitude               = 60.0,  # degrees latitude
    equilibrium_equatorial_temperature           = 25.0,  # degrees Celsius
    equilibrium_polar_temperature                = 2.0,  # degrees Celsius
    initial_temperature_scale_height             = 400.0,
    initial_bottom_temperature                   = 0.0,
    initial_velocity_noise_amplitude             = 1e-2,
    bottom_drag_coefficient                      = 3e-3,
    equation_of_state                            = LinearEquationOfState(thermal_expansion=2e-4, haline_contraction=8e-5),
    # equation_of_state                            = TEOS10EquationOfState(; reference_density)
    time_step                                    = 20minutes,
    stop_iteration                               = Inf,
    stop_time                                    = Inf,
    tracers                                      = [:T, :S, :e])

    # A spherical domain
    grid = LatitudeLongitudeGrid(architecture;
                                 size, longitude, latitude, z,
                                 halo = (5, 5, 5))
    
    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)

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

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)

    planet_extent = latitude[2] - latitude[1]

    wind_stress_parameters = (; maximum_wind_stress,
                              equatorial_wind_extent,
                              planet_extent)

    u_wind_stress_bc = FluxBoundaryCondition(aqua_planet_wind_stress, discrete_form=true, parameters=wind_stress_parameters)

    T_relaxation_parameters = (; surface_temperature_relaxation_time_scale,
                                 temperature_asymmetry_factor,
                                 temperature_asymmetry_latitude,
                                 equilibrium_equatorial_temperature,
                                 equilibrium_polar_temperature,
                                 planet_extent)
     
    T_surface_relaxation_bc = FluxBoundaryCondition(aqua_planet_temperature_flux,
                                                    discrete_form = true,
                                                    parameters = T_relaxation_parameters)

    u_bcs = FieldBoundaryConditions(top = u_wind_stress_bc,
                                    bottom = u_bottom_drag_bc,
                                    immersed = u_immersed_bc)

    v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc,
                                    immersed = v_immersed_bc)

    T_bcs = FieldBoundaryConditions(top = T_surface_relaxation_bc)

    buoyancy = SeawaterBuoyancy(; equation_of_state)
    coriolis = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme())
    free_surface = ImplicitFreeSurface()

    @info "Building a model..."; start=time_ns()

    model = HydrostaticFreeSurfaceModel(; grid, free_surface, buoyancy, coriolis, tracers, closure,
                                        momentum_advection = VectorInvariant(vorticity_scheme   = WENO(),
                                                                             divergence_scheme  = WENO(),
                                                                             vertical_scheme    = WENO()),
                                        tracer_advection = WENO(grid),
                                        boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs))

    @info "... built $model."
    @info "Model building time: " * prettytime(1e-9 * (time_ns() - start))

    #####
    ##### Initial condition
    #####

    Tᵢ = AquaPlanetInitialCondition(initial_temperature_scale_height,
                                    initial_bottom_temperature,
                                    T_relaxation_parameters)

    ℓ = 10 # degrees
    ϵ = 100.0 # m s⁻¹
    k = 8 * π/planet_extent # degrees
    uᵢ(x, y, z) = ϵ * ũ(x, y, ℓ, k) + initial_velocity_noise_amplitude * (2rand() - 1)
    vᵢ(x, y, z) = ϵ * ṽ(x, y, ℓ, k) + initial_velocity_noise_amplitude * (2rand() - 1)

    set!(model, u=uᵢ, v=vᵢ, T=Tᵢ, e=1e-6, S=35.0)

    simulation = Simulation(model; Δt=time_step, stop_iteration, stop_time)

    start_time = Ref(time_ns())

    function progress(sim)
        wall_time = (time_ns() - start_time[]) * 1e-9

        u = sim.model.velocities.u
        v = sim.model.velocities.v
        w = sim.model.velocities.w

        T = sim.model.tracers.T
        S = sim.model.tracers.S
        e = sim.model.tracers.e

        u_interior = Array(interior(u))
        v_interior = Array(interior(v))
        w_interior = Array(interior(w))
        T_interior = Array(interior(T))
        S_interior = Array(interior(S))
        e_interior = Array(interior(e))

        max_w, i_max_w = findmax(w_interior)
        max_u, i_max_u = findmax(u_interior)
        max_v, i_max_v = findmax(v_interior)
        max_T, i_max_T = findmax(T_interior)
        max_S, i_max_S = findmax(S_interior)
        max_e, i_max_e = findmax(e_interior)

        msg1 = @sprintf("Time: % 12s, iteration: %d, ", prettytime(sim), iteration(sim))

        msg2 = @sprintf("max(u): %.2e (%d, %d, %d) m s⁻¹, max(v): %.2e (%d, %d, %d) m s⁻¹, max(w): %.2e (%d, %d, %d) m s⁻¹, ",
                        max_u, i_max_u[1], i_max_u[2], i_max_u[3],
                        max_v, i_max_v[1], i_max_v[2], i_max_v[3],
                        max_w, i_max_w[1], i_max_w[2], i_max_w[3])

        msg3 = @sprintf("max(T): %.2e (%d, %d, %d) ᵒC, ",
                        max_T, i_max_T[1], i_max_T[2], i_max_T[3])

        msg4 = @sprintf("extrema(e): (%.2e, %.2e)  m² s⁻², (%d, %d, %d), ",
                        maximum(e), minimum(e), i_max_e[1], i_max_e[2], i_max_e[3])

        msg5 = @sprintf("wall time: %s", prettytime(wall_time))

        @info msg1 * msg2 * msg3 * msg4 * msg5

        start_time[] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    return simulation
end

