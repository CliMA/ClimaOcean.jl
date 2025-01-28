using Oceananigans.BuoyancyFormulations: g_Earth

struct CoefficientBasedFluxes{CD, CH, CQ, ΔU, FT}
    drag_coefficient :: CD
    gravitational_acceleration :: FT
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
                                gravitational_acceleration = g_Earth,
                                heat_transfer_coefficient = drag_coefficient,
                                vapor_flux_coefficient = drag_coefficient,
                                bulk_velocity = RelativeVelocity(),
                                solver_tolerance = 1e-8,
                                solver_maxiter = 100)

    drag_coefficient = convert_if_number(FT, drag_coefficient)
    heat_transfer_coefficient = convert_if_number(FT, heat_transfer_coefficient)
    vapor_flux_coefficient = convert_if_number(FT, vapor_flux_coefficient)

    return CoefficientBasedFluxes(drag_coefficient,
                                  gravitational_acceleration,
                                  heat_transfer_coefficient,
                                  vapor_flux_coefficient,
                                  bulk_velocity,
                                  convert(FT, solver_tolerance),
                                  solver_maxiter)
end

function iterate_interface_fluxes(flux_formulation::CoefficientBasedFluxes,
                                  Tₛ, qₛ, Δθ, Δq, Δh,
                                  approximate_interface_state,
                                  atmosphere_state,
                                  atmosphere_properties)

    Ψₐ = atmosphere_state
    Ψ̃ᵢ = approximate_interface_state
    Δu, Δv = velocity_difference(flux_formulation.bulk_velocity, Ψₐ, Ψ̃ᵢ)
    ΔU = sqrt(Δu^2 + Δv^2)

    Cd = flux_formulation.drag_coefficient
    Ch = flux_formulation.heat_transfer_coefficient
    Cq = flux_formulation.vapor_flux_coefficient

    u★ = sqrt(Cd) * ΔU
    # u★ θ★ = Ch * Δθ * ΔU
    θ★ = Ch / sqrt(Cd) * Δθ
    q★ = Cq / sqrt(Cd) * Δq

    return u★, θ★, q★
end

