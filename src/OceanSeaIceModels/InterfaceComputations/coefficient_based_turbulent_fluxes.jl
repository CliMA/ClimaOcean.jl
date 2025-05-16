using Oceananigans.BuoyancyFormulations: g_Earth

"""
    CoefficientBasedFluxes(FT = Oceananigans.defaults.FloatType;
                          drag_coefficient = 1e-3,
                          heat_transfer_coefficient = drag_coefficient,
                          vapor_flux_coefficient = drag_coefficient,
                          solver_stop_criteria = nothing,
                          solver_tolerance = 1e-8,
                          solver_maxiter = 20)

A structure for computing turbulent fluxes using constant bulk transfer coefficients.

# Arguments
- `FT`: (optional) Float type for the coefficients, defaults to Oceananigans.defaults.FloatType
- `drag_coefficient`: Coefficient for momentum transfer (Cᵈ), defaults to 1e-3
- `heat_transfer_coefficient`: Coefficient for heat transfer (Cʰ), defaults to drag_coefficient
- `vapor_flux_coefficient`: Coefficient for moisture transfer (Cᵍ), defaults to drag_coefficient
- `solver_stop_criteria`: Criteria for iterative solver convergence. If nothing, creates new criteria using tolerance and maxiter
- `solver_tolerance`: Tolerance for solver convergence when creating new stop criteria, defaults to 1e-8
- `solver_maxiter`: Maximum iterations for solver when creating new stop criteria, defaults to 20

# Fields
- `drag_coefficient`: Coefficient for momentum transfer
- `heat_transfer_coefficient`: Coefficient for sensible heat transfer
- `vapor_flux_coefficient`: Coefficient for latent heat transfer
- `solver_stop_criteria`: Criteria for iterative solver convergence

Used in bulk flux calculations to determine the exchange of momentum, heat, and moisture 
between the ocean/ice surface and the atmosphere using constant transfer coefficients.

# Example
```jldoctest
using Oceananigans
using ClimaOcean

grid = RectilinearGrid(size=3, z=(-1, 0), topology=(Flat, Flat, Bounded))
ocean = ocean_simulation(grid)

ao_fluxes = CoefficientBasedFluxes(drag_coefficient = 1e-2,
                                   heat_transfer_coefficient = 1e-3,
                                   vapor_flux_coefficient = 1e-3)

interfaces = ComponentInterfaces(nothing, ocean; atmosphere_ocean_flux_formulation=ao_fluxes)

# output
ComponentInterfaces
```
"""
struct CoefficientBasedFluxes{CD, CH, CQ, S}
    drag_coefficient :: CD
    heat_transfer_coefficient :: CH
    vapor_flux_coefficient :: CQ
    solver_stop_criteria :: S
end

convert_if_number(FT, a::Number) = convert(FT, a)
convert_if_number(FT, a) = a

function CoefficientBasedFluxes(FT = Oceananigans.defaults.FloatType;
                                drag_coefficient = 1e-3,
                                heat_transfer_coefficient = drag_coefficient,
                                vapor_flux_coefficient = drag_coefficient,
                                solver_stop_criteria = nothing,
                                solver_tolerance = 1e-8,
                                solver_maxiter = 20)

    if isnothing(solver_stop_criteria)
        solver_tolerance = convert(FT, solver_tolerance)
        solver_stop_criteria = ConvergenceStopCriteria(solver_tolerance, solver_maxiter)
    end

    drag_coefficient = convert_if_number(FT, drag_coefficient)
    heat_transfer_coefficient = convert_if_number(FT, heat_transfer_coefficient)
    vapor_flux_coefficient = convert_if_number(FT, vapor_flux_coefficient)

    return CoefficientBasedFluxes(drag_coefficient,
                                  heat_transfer_coefficient,
                                  vapor_flux_coefficient,
                                  solver_stop_criteria)
end

@inline function iterate_interface_fluxes(flux_formulation::CoefficientBasedFluxes,
                                          Tₛ, qₛ, Δθ, Δq, Δh,
                                          approximate_interface_state,
                                          atmosphere_state,
                                          interface_properties,
                                          atmosphere_properties)

    Ψₐ = atmosphere_state
    Ψ̃ᵢ = approximate_interface_state
    Δu, Δv = velocity_difference(interface_properties.velocity_formulation, Ψₐ, Ψ̃ᵢ)
    ΔU = sqrt(Δu^2 + Δv^2)

    Cd = flux_formulation.drag_coefficient
    Ch = flux_formulation.heat_transfer_coefficient
    Cq = flux_formulation.vapor_flux_coefficient

    u★ = sqrt(Cd) * ΔU
    θ★ = Ch / sqrt(Cd) * Δθ
    q★ = Cq / sqrt(Cd) * Δq

    return u★, θ★, q★
end

