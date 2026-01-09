#####
##### Iterating conditions for stop criteria
#####

@inline function iterating(Ψⁿ, Ψ⁻, iteration, convergence::ConvergenceStopCriteria)
    maxiter = convergence.maxiter
    tolerance = convergence.tolerance
    hasnt_started = iteration == 0
    reached_maxiter = iteration ≥ maxiter
    drift = abs(Ψⁿ.u★ - Ψ⁻.u★) + abs(Ψⁿ.θ★ - Ψ⁻.θ★) + abs(Ψⁿ.q★ - Ψ⁻.q★)
    converged = drift < tolerance
    return !(converged | reached_maxiter) | hasnt_started
end

@inline iterating(Ψⁿ, Ψ⁻, iteration, fixed::FixedIterations) = iteration < fixed.iterations

#####
##### Main solver dispatch
#####

"""
    compute_interface_state(flux_formulation, initial_interface_state, ...)

Compute the interface state (u★, θ★, q★, T, ...) by iteratively solving
the similarity theory equations. Dispatches to the appropriate solver
based on `flux_formulation.solver`.

Returns the converged `InterfaceState`.
"""
@inline function compute_interface_state(flux_formulation,
                                         initial_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties)

    solver = flux_formulation.solver
    return compute_interface_state(solver,
                                   flux_formulation,
                                   initial_interface_state,
                                   atmosphere_state,
                                   interior_state,
                                   downwelling_radiation,
                                   interface_properties,
                                   atmosphere_properties,
                                   interior_properties)
end

#####
##### Fixed-point iteration solver
#####

"""
    compute_interface_state(solver::FixedPointSolver, flux_formulation, ...)

Compute the interface state using fixed-point (Picard) iteration.

The iteration proceeds by repeatedly applying `iterate_interface_state`
until the characteristic scales (u★, θ★, q★) converge or the maximum
number of iterations is reached.
"""
@inline function compute_interface_state(solver::FixedPointSolver,
                                         flux_formulation,
                                         initial_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties)

    Ψₐ = atmosphere_state
    Ψᵢ = interior_state
    Ψₛⁿ = Ψₛ⁻ = initial_interface_state
    stop_criteria = solver.stop_criteria
    iteration = 0

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

    return Ψₛⁿ
end

#####
##### Broyden's method solver
#####

"""
    compute_interface_state(solver::BroydenSolver{N}, flux_formulation, ...)

Compute the interface state using Good Broyden's quasi-Newton method.

This method typically converges faster than fixed-point iteration by
approximating the Jacobian of the residual function and using rank-1
updates (Sherman-Morrison formula).

For N=3 (BulkTemperature): solves for (u★, θ★, q★)
For N=4 (SkinTemperature): solves for (u★, θ★, q★, Tₛ)
"""
@inline function compute_interface_state(solver::BroydenSolver{N, FT},
                                         flux_formulation,
                                         initial_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties) where {N, FT}

    Ψₐ = atmosphere_state
    Ψᵢ = interior_state
    Ψₛ = initial_interface_state

    stop_criteria = solver.stop_criteria
    tolerance = stop_criteria.tolerance
    maxiter = stop_criteria.maxiter

    valN = Val(N)

    # Initialize inverse Jacobian as scaled identity
    J_inv = identity_matrix_inv(valN, solver.initial_jacobian_scale)

    # Compute initial residual F = iterate(x) - x
    F, Ψₛ_new = compute_broyden_residual(flux_formulation,
                                         Ψₛ, Ψₐ, Ψᵢ,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties,
                                         valN)

    x = state_to_tuple(Ψₛ, valN)
    iteration = 0
    converged = broyden_converged(F, tolerance, valN)

    # Main Broyden iteration loop
    # Using bitwise & instead of && for GPU compatibility
    while !converged & (iteration < maxiter)
        # Compute Newton-like step: Δx = -J⁻¹ F
        Δx = negative_mat_vec_mul(J_inv, F, valN)

        # Update state: x_new = x + Δx
        x_new = add_tuples(x, Δx, valN)

        # Create new interface state from updated tuple
        Ψₛ = tuple_to_state(x_new, Ψₛ, valN)

        # Store old residual
        F_old = F

        # Compute new residual
        F, Ψₛ_new = compute_broyden_residual(flux_formulation,
                                             Ψₛ, Ψₐ, Ψᵢ,
                                             downwelling_radiation,
                                             interface_properties,
                                             atmosphere_properties,
                                             interior_properties,
                                             valN)

        # Compute residual change
        ΔF = subtract_tuples(F, F_old, valN)

        # Update inverse Jacobian using Good Broyden formula
        J_inv = broyden_update(J_inv, Δx, ΔF, valN)

        # Update for next iteration
        x = x_new
        iteration += 1
        converged = broyden_converged(F, tolerance, valN)
    end

    # Return the final state after one more iteration for consistency
    # This ensures the returned state has properly computed T, q, etc.
    return Ψₛ_new
end

#####
##### Single iteration step
#####

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

    # Recompute the saturation specific humidity at the interface based on the new temperature
    q_formulation = interface_properties.specific_humidity_formulation
    Sₛ = approximate_interface_state.S
    qₛ = surface_specific_humidity(q_formulation, ℂₐ, 𝒬ₐ, Tₛ, Sₛ)

    # Compute the specific humidity increment
    qₐ = AtmosphericThermodynamics.vapor_specific_humidity(ℂₐ, 𝒬ₐ)
    Δq = qₐ - qₛ

    θₐ = surface_atmosphere_temperature(atmosphere_state, atmosphere_properties)
    Δθ = θₐ - Tₛ
    Δh = atmosphere_state.z # Assumption! The surface is at z = 0 -> Δh = zₐ - 0

    u★, θ★, q★ = iterate_interface_fluxes(flux_formulation,
                                          Tₛ, qₛ, Δθ, Δq, Δh,
                                          approximate_interface_state,
                                          atmosphere_state,
                                          interface_properties,
                                          atmosphere_properties)

    u = approximate_interface_state.u
    v = approximate_interface_state.v
    S = approximate_interface_state.S

    return InterfaceState(convert(FT, u★),
                          convert(FT, θ★),
                          convert(FT, q★),
                          convert(FT, u),
                          convert(FT, v),
                          convert(FT, Tₛ),
                          convert(FT, S),
                          convert(FT, qₛ))
end
