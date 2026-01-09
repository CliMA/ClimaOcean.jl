using DocStringExtensions

"""
    struct CoefficientBasedFluxes{CD, CH, CQ, S}

A structure for computing turbulent fluxes using constant bulk transfer coefficients.

$(TYPEDFIELDS)
"""
struct CoefficientBasedFluxes{CD, CH, CQ, S}
    "Coefficient for momentum transfer"
    drag_coefficient :: CD
    "Coefficient for sensible heat transfer"
    heat_transfer_coefficient :: CH
    "Coefficient for latent heat transfer"
    vapor_flux_coefficient :: CQ
    "Iterative solver for computing interface state"
    solver :: S
end

Base.summary(flux_formulation::CoefficientBasedFluxes) = "CoefficientBasedFluxes"

function Base.show(io::IO, flux_formulation::CoefficientBasedFluxes)
    print(io, summary(flux_formulation), '\n')
    print(io, "├── drag_coefficient: ", prettysummary(flux_formulation.drag_coefficient), '\n')
    print(io, "├── heat_transfer_coefficient: ", prettysummary(flux_formulation.heat_transfer_coefficient), '\n')
    print(io, "├── vapor_flux_coefficient: ", prettysummary(flux_formulation.vapor_flux_coefficient), '\n')
    print(io, "└── solver: ", summary(flux_formulation.solver), '\n')
end

convert_if_number(FT, a::Number) = convert(FT, a)
convert_if_number(FT, a) = a

"""
    CoefficientBasedFluxes(FT = Oceananigans.defaults.FloatType;
                          drag_coefficient = 1e-3,
                          heat_transfer_coefficient = drag_coefficient,
                          vapor_flux_coefficient = drag_coefficient,
                          solver = nothing,
                          solver_tolerance = 1e-8,
                          solver_maxiter = 20)

Return the structure for computing turbulent fluxes using constant bulk transfer coefficients.
Used in bulk flux calculations to determine the exchange of momentum, heat, and moisture
between the ocean/ice surface and the atmosphere using constant transfer coefficients.

# Arguments
- `FT`: (optional) Float type for the coefficients, defaults to Oceananigans.defaults.FloatType

# Keyword Arguments
- `drag_coefficient`: Coefficient for momentum transfer (`Cᵈ`), defaults to 1e-3.
- `heat_transfer_coefficient`: Coefficient for heat transfer (`Cʰ`), defaults to `drag_coefficient`
- `vapor_flux_coefficient`: Coefficient for moisture transfer (`Cᵍ`), defaults to `drag_coefficient`
- `solver`: The iterative solver to use. Default: `FixedPointSolver`.
            Can also use `BroydenSolver(Val(N))` for quasi-Newton acceleration.
- `solver_tolerance`: Tolerance for solver convergence (used if `solver` is not specified), defaults to 1e-8.
- `solver_maxiter`: Maximum iterations for solver (used if `solver` is not specified), defaults to 20

# Example
```jldoctest
using Oceananigans
using ClimaOcean

grid = RectilinearGrid(size=3, z=(-1, 0), topology=(Flat, Flat, Bounded))
ocean = ocean_simulation(grid; timestepper = :QuasiAdamsBashforth2)

ao_fluxes = CoefficientBasedFluxes(drag_coefficient = 1e-2,
                                   heat_transfer_coefficient = 1e-3,
                                   vapor_flux_coefficient = 1e-3)

interfaces = ComponentInterfaces(nothing, ocean; atmosphere_ocean_fluxes=ao_fluxes)

# output
ComponentInterfaces
```
"""
function CoefficientBasedFluxes(FT = Oceananigans.defaults.FloatType;
                                drag_coefficient = 1e-3,
                                heat_transfer_coefficient = drag_coefficient,
                                vapor_flux_coefficient = drag_coefficient,
                                solver = nothing,
                                solver_stop_criteria = nothing,  # Deprecated, for backwards compatibility
                                solver_tolerance = 1e-8,
                                solver_maxiter = 20)

    # Handle solver construction with backwards compatibility
    if isnothing(solver)
        if isnothing(solver_stop_criteria)
            solver_tolerance = convert(FT, solver_tolerance)
            solver_stop_criteria = ConvergenceStopCriteria(solver_tolerance, solver_maxiter)
        end
        solver = FixedPointSolver(solver_stop_criteria)
    end

    drag_coefficient = convert_if_number(FT, drag_coefficient)
    heat_transfer_coefficient = convert_if_number(FT, heat_transfer_coefficient)
    vapor_flux_coefficient = convert_if_number(FT, vapor_flux_coefficient)

    return CoefficientBasedFluxes(drag_coefficient,
                                  heat_transfer_coefficient,
                                  vapor_flux_coefficient,
                                  solver)
end

# Backwards compatibility: allow access to solver_stop_criteria
Base.getproperty(f::CoefficientBasedFluxes, name::Symbol) =
    name === :solver_stop_criteria ? getfield(f, :solver).stop_criteria : getfield(f, name)

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
