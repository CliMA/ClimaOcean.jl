

struct CoefficientBasedFluxes{CD, CH, CQ, ΔU, FT}
    drag_coefficient :: CD
    heat_transfer_coefficient :: CH
    vapor_flux_coefficient :: CQ
    bulk_velocity :: ΔU
    solver_tolerance :: FT
    solver_maxiter :: Int
end

convert_if_number(FT, a::Number) = convert(FT, a)
convert_if_number(FT, a) = a

function CoefficientBasedFluxes(FT = Float64;
                                drag_coefficient = 1e-3,
                                heat_transfer_coefficient = drag_coefficient,
                                vapor_flux_coefficient = drag_coefficient,
                                bulk_velocity = RelativeVelocity(),
                                solver_tolerance = 1e-8,
                                solver_maxiter = 10)

    drag_coefficient = convert_if_number(FT, drag_coefficient)
    heat_transfer_coefficient = convert_if_number(FT, heat_transfer_coefficient)
    vapor_flux_coefficient = convert_if_number(FT, vapor_flux_coefficient)

    return CoefficientBasedFluxes(drag_coefficient,
                                  heat_transfer_coefficient,
                                  vapor_flux_coefficient,
                                  bulk_velocity,
                                  convert(FT, solver_tolerance),
                                  solver_maxiter)
end

@inline function compute_interface_state(flux_formulation::CoefficientBasedFluxes,
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
    iteration = 0
    maxiter = flux_formulation.solver_maxiter
    tolerance = flux_formulation.solver_tolerance

    while iterating(Ψₛⁿ, Ψₛ⁻, iteration, maxiter, tolerance)
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

"""
    iterate_interface_state(flux_formulation, Ψₛⁿ⁻¹, Ψₐ, Ψᵢ, Qᵣ, ℙₛ, ℙₐ, ℙᵢ)

Return the nth iterate of the interface state `Ψₛⁿ` computed according to the
`flux_formulation`, given the interface state at the previous iterate `Ψₛⁿ⁻¹`,
as well as the atmosphere state `Ψₐ`, the interior state `Ψᵢ`,
downwelling radiation `Qᵣ`, and the interface, atmosphere,
and interior properties `ℙₛ`, `ℙₐ`, and `ℙᵢ`.
"""
@inline function iterate_interface_state(flux_formulation::CoefficientBasedFluxes,
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
    g = 9.81 #flux_formulation.gravitational_acceleration
    cₐ = interior_properties.heat_capacity
    θₐ = Tₐ + g * Δh / cₐ
    Δθ = θₐ - Tₛ

    Δu, Δv = velocity_difference(flux_formulation.bulk_velocity, atmosphere_state, approximate_interface_state)
    ΔU = sqrt(Δu^2 + Δv^2)

    Cd = flux_formulation.drag_coefficient
    Ch = flux_formulation.heat_transfer_coefficient
    Cq = flux_formulation.vapor_flux_coefficient

    u★ = sqrt(Cd) * ΔU
    θ★ = Ch / sqrt(Cd) * Δθ
    q★ = Cq / sqrt(Cd) * Δq

    u = approximate_interface_state.u
    v = approximate_interface_state.v
    S = approximate_interface_state.S

    return InterfaceState(u★, θ★, q★, u, v, Tₛ, S, convert(FT, qₛ))
end

