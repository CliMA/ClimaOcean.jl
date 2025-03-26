#####
##### Solver stop criteria
#####

abstract type AbstractIterativeSolver{S} end

struct ConvergenceStopCriteria{FT}
    tolerance :: FT
    maxiter :: Int
end

@inline function iterating(Ψⁿ, Ψ⁻, iteration, convergence::ConvergenceStopCriteria)
    maxiter = convergence.maxiter
    tolerance = convergence.tolerance
    hasnt_started = iteration == 0
    reached_maxiter = iteration ≥ maxiter
    drift = abs(Ψⁿ.u★ - Ψ⁻.u★) + abs(Ψⁿ.θ★ - Ψ⁻.θ★) + abs(Ψⁿ.q★ - Ψ⁻.q★)
    converged = drift < tolerance
    return !(converged | reached_maxiter) | hasnt_started
end

struct FixedIterations{I}
    iterations :: I
end

@inline iterating(Ψⁿ, Ψ⁻, iteration, fixed::FixedIterations) = iteration < fixed.iterations

#####
##### The solver
#####

# Iterating condition for the characteristic scales solvers
@inline function compute_interface_state(flux_formulation::AbstractIterativeSolver,
                                         need_to_iterate, # Bool -> do we need to iterate or not?
                                         initial_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties)


    stop_criteria = flux_formulation.solver_stop_criteria
    needs_to_converge = stop_criteria isa ConvergenceStopCriteria

    Ψₐ = atmosphere_state
    Ψᵢ = interior_state
    Ψₛⁿ = Ψₛ⁻ = initial_interface_state
    stop_criteria = flux_formulation.solver_stop_criteria
    iteration = 0

    if needs_to_converge && !need_to_iterate
        Ψₛⁿ = zero_interface_state(FT)
    else
        while iterating(Ψₛⁿ, Ψₛ⁻, iteration, stop_criteria)
            Ψₛ⁻ = Ψₛⁿ
            Ψₛⁿ = iterate_interface_state(flux_formulation,
                                          Ψₛ⁻, Ψₐ, Ψᵢ,
                                          downwelling_radiation,
                                          interface_properties,
                                          atmosphere_properties,
                                          interior_properties)
            iteration += 1
        end
    end

    return Ψₛⁿ
end

"""
    iterate_interface_state(flux_formulation, Ψₛⁿ⁻¹, Ψₐ, Ψᵢ, Qᵣ, ℙₛ, ℙₐ, ℙᵢ)

Return the nth iterate of the interface state `Ψₛⁿ` computed according to the
`flux_formulation`, given the interface state at the previous iterate `Ψₛⁿ⁻¹`,
as well as the atmosphere state `Ψₐ`, the interior state `Ψᵢ`,
downwelling radiation `Qᵣ`, and the interface, atmosphere,
and interior properties `ℙₛ`, `ℙₐ`, and `ℙᵢ`.
"""
@inline function iterate_interface_state(flux_formulation,
                                         approximate_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties)

    Tₛ = compute_interface_temperature(interface_properties.temperature_formulation,
                                       approximate_interface_state,
                                       atmosphere_state,
                                       interior_state,
                                       downwelling_radiation,
                                       interface_properties,
                                       atmosphere_properties,
                                       interior_properties)

    # Thermodynamic state
    FT = eltype(approximate_interface_state)
    ℂₐ = atmosphere_properties.thermodynamics_parameters
    𝒬ₐ = atmosphere_state.𝒬
    ρₐ = 𝒬ₐ.ρ

    # Recompute the saturation specific humidity at the interface based on the new temperature
    q_formulation = interface_properties.specific_humidity_formulation
    Sₛ = approximate_interface_state.S
    qₛ = saturation_specific_humidity(q_formulation, ℂₐ, ρₐ, Tₛ, Sₛ)

    # Compute the specific humidity increment
    qₐ = AtmosphericThermodynamics.vapor_specific_humidity(ℂₐ, 𝒬ₐ)
    Δq = qₐ - qₛ

    # Temperature increment including the ``lapse rate'' `α = g / cₚ`
    zₐ = atmosphere_state.z
    zₛ = zero(FT)
    Δh = zₐ - zₛ
    Tₐ = AtmosphericThermodynamics.air_temperature(ℂₐ, 𝒬ₐ)
    g  = flux_formulation.gravitational_acceleration
    cₐ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ)
    θₐ = Tₐ + g * Δh / cₐ
    Δθ = θₐ - Tₛ

    u★, θ★, q★ = iterate_interface_fluxes(flux_formulation,
                                          Tₛ, qₛ, Δθ, Δq, Δh,
                                          approximate_interface_state,
                                          atmosphere_state,
                                          interface_properties,
                                          atmosphere_properties)

    u = approximate_interface_state.u
    v = approximate_interface_state.v
    S = approximate_interface_state.S

    return InterfaceState(u★, θ★, q★, u, v, Tₛ, S, convert(FT, qₛ))
end
