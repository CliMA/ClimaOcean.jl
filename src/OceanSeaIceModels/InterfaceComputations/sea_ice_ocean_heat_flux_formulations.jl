using ClimaSeaIce.SeaIceThermodynamics: melting_temperature, LinearLiquidus

#####
##### Ice Bath Heat Flux (bulk formulation)
#####

"""
    IceBathHeatFlux{FT, U}

Bulk formulation for sea ice-ocean heat flux.

The interface temperature is fixed at the freezing point of the surface salinity,
and the heat flux is computed using bulk transfer:
```math
Q = \\rho_o c_o \\alpha_h u_* (T - T_m)
```
where ``\\alpha_h`` is the heat transfer coefficient and ``u_*`` is the friction velocity.

Fields
======

- `heat_transfer_coefficient::FT`: turbulent heat exchange coefficient ``\\alpha_h`` (dimensionless)
- `friction_velocity::U`: friction velocity value or formulation (constant `Number` or `MomentumBasedFrictionVelocity`)

Example
=======

```jldoctest
using ClimaOcean.OceanSeaIceModels: IceBathHeatFlux

flux = IceBathHeatFlux(heat_transfer_coefficient = 0.006, friction_velocity = 0.002)

# output
IceBathHeatFlux{Float64}
├── heat_transfer_coefficient: 0.006
└── friction_velocity: 0.002
```

References
==========

- [holland1999modeling](@citet): Holland, D. M., & Jenkins, A. (1999). Modeling thermodynamic ice–ocean interactions
  at the base of an ice shelf. *Journal of Physical Oceanography*, 29(8), 1787-1800.
"""
struct IceBathHeatFlux{FT, U}
    heat_transfer_coefficient :: FT
    friction_velocity :: U
end

"""
    IceBathHeatFlux(FT::DataType = Float64;
                    heat_transfer_coefficient = 0.006,
                    friction_velocity = 0.02)

Construct an `IceBathHeatFlux` with the specified parameters.

Keyword Arguments
=================

- `heat_transfer_coefficient`: turbulent heat exchange coefficient. Default: 0.006.
- `friction_velocity`: friction velocity value or formulation. Default: 0.02.
"""
function IceBathHeatFlux(FT::DataType = Float64;
                         heat_transfer_coefficient = convert(FT, 0.006),
                         friction_velocity = convert(FT, 0.02))
    return IceBathHeatFlux(convert(FT, heat_transfer_coefficient), friction_velocity)
end

#####
##### Three-Equation Heat Flux (full formulation)
#####

"""
    ThreeEquationHeatFlux{FT, U}

Three-equation formulation for sea ice-ocean heat flux.

This formulation solves a coupled system for the interface temperature and salinity:
1. Heat balance: ``\\rho c_p \\gamma_T (T - T_b) = ℰ q``
2. Salt balance: ``\\gamma_S (S - S_b) = q (S_b - S_i)``
3. Freezing point: ``T_b = T_m(S_b)``

where ``T_b`` and ``S_b`` are the interface temperature and salinity,
``\\gamma_T = \\alpha_h u_*`` and ``\\gamma_S = \\alpha_s u_*`` are turbulent exchange velocities,
``L`` is the latent heat of fusion, and ``q`` is the melt rate (computed, not input).

Fields
======

- `heat_transfer_coefficient::FT`: turbulent heat exchange coefficient ``\\alpha_h`` (dimensionless)
- `salt_transfer_coefficient::FT`: turbulent salt exchange coefficient ``\\alpha_s`` (dimensionless)
- `friction_velocity::U`: friction velocity value or formulation (constant `Number` or `MomentumBasedFrictionVelocity`)

Example
=======

```jldoctest
using ClimaOcean.OceanSeaIceModels: ThreeEquationHeatFlux

flux = ThreeEquationHeatFlux()

# output
ThreeEquationHeatFlux{Float64}
├── heat_transfer_coefficient: 0.0095
├── salt_transfer_coefficient: 0.00027142857142857146
└── friction_velocity: 0.002
```

References
==========

- [holland1999modeling](@citet): Holland, D. M., & Jenkins, A. (1999). Modeling thermodynamic ice–ocean interactions
  at the base of an ice shelf. *Journal of Physical Oceanography*, 29(8), 1787-1800.
- [hieronymus2021comparison](@citet): Hieronymus, M., et al. (2021). A comparison of ocean-ice flux parametrizations.
  *Geosci. Model Dev.*, 14, 4891-4908.
"""
struct ThreeEquationHeatFlux{FT, U}
    heat_transfer_coefficient :: FT
    salt_transfer_coefficient :: FT
    friction_velocity :: U
end

"""
    ThreeEquationHeatFlux(FT::DataType = Float64;
                          heat_transfer_coefficient = 0.0095,
                          salt_transfer_coefficient = heat_transfer_coefficient / 35,
                          friction_velocity = 0.002)

Construct a `ThreeEquationHeatFlux` with the specified parameters.

Default values follow [hieronymus2021comparison](@citet) with ``R = \\alpha_h / \\alpha_s = 35``.

Keyword Arguments
=================

- `heat_transfer_coefficient`: turbulent heat exchange coefficient ``\\alpha_h``. Default: 0.0095.
- `salt_transfer_coefficient`: turbulent salt exchange coefficient ``\\alpha_s``. Default: ``\\alpha_h / 35 \\approx 0.000271``.
- `friction_velocity`: friction velocity value or formulation. Default: 0.002.
"""
function ThreeEquationHeatFlux(FT::DataType = Oceananigans.defaults.FloatType;
                               heat_transfer_coefficient = 0.0095,
                               salt_transfer_coefficient = heat_transfer_coefficient / 35,
                               friction_velocity = convert(FT, 0.002))
    return ThreeEquationHeatFlux(convert(FT, heat_transfer_coefficient),
                                 convert(FT, salt_transfer_coefficient),
                                 friction_velocity)
end

#####
##### Interface heat flux computation
#####

"""
    compute_interface_heat_flux(flux::IceBathHeatFlux, ...)

Compute the heat flux and melt rate at the sea ice-ocean interface using bulk formulation.
Returns `(Q, q)` where:
- `Q > 0` means heat flux from ocean to ice (ocean cooling)
- `q > 0` means melting (ice volume loss)
"""
@inline function compute_interface_heat_flux(flux::IceBathHeatFlux, i, j,
                                              Tᵢ, Sᵢ, Tₒ, Sₒ, Sⁱ, ℵ, Nz,
                                              liquidus, ρₒ, cₒ, ℰ, τx, τy)
    @inbounds begin
        Tᴺ  = Tᵢ[i, j, 1]
        Sᴺ  = Sᵢ[i, j, 1]
        ℵᵢⱼ = ℵ[i, j, 1]
    end

    # Interface temperature is at the freezing point
    Tₘ = melting_temperature(liquidus, Sᴺ)

    αₕ = flux.heat_transfer_coefficient
    u★ = get_friction_velocity(flux.friction_velocity, i, j, τx, τy, ρₒ)

    # Heat flux: Q > 0 means heat flux from ocean to ice (ocean cooling)
    Qᵢₒ = ρₒ * cₒ * αₕ * u★ * (Tᴺ - Tₘ) * ℵᵢⱼ

    # Melt rate: q = Q / L (positive for melting)
    q = Qᵢₒ / ℰ

    return Qᵢₒ, q
end

"""
    compute_interface_heat_flux(flux::ThreeEquationHeatFlux, ...)

Compute the heat flux and melt rate at the sea ice-ocean interface using three-equation formulation.
Returns `(Q, q)` where:
- `Q > 0` means heat flux from ocean to ice (ocean cooling)
- `q > 0` means melting (ice volume loss)
"""
@inline function compute_interface_heat_flux(flux::ThreeEquationHeatFlux, i, j,
                                              Tᵢ, Sᵢ, Tₒ, Sₒ, Sⁱ, ℵ, Nz,
                                              liquidus, ρₒ, cₒ, ℰ, τx, τy)
    @inbounds begin
        Tᴺ   = Tₒ[i, j, Nz]   # Ocean surface temperature
        Sᴺ   = Sₒ[i, j, Nz]   # Ocean surface salinity
        Sᵢᶜᵉ = Sⁱ[i, j, 1]    # Ice salinity
        ℵᵢⱼ  = ℵ[i, j, 1]     # Ice concentration
    end

    αₕ = flux.heat_transfer_coefficient
    αₛ = flux.salt_transfer_coefficient
    u★ = get_friction_velocity(flux.friction_velocity, i, j, τx, τy, ρₒ)

    # Solve for interface temperature, salinity, and melt rate
    Tᵦ, Sᵦ, q = solve_interface_conditions(Tᴺ, Sᴺ, Sᵢᶜᵉ, αₕ, αₛ, u★, ℰ, ρₒ, cₒ, liquidus)

    # Store interface values
    @inbounds Tᵢ[i, j, 1] = Tᵦ
    @inbounds Sᵢ[i, j, 1] = Sᵦ

    # Heat flux: Q > 0 means heat flux from ocean to ice (ocean cooling)
    Qᵢₒ = ℰ * q * ℵᵢⱼ

    return Qᵢₒ, q
end

"""
    solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★, ℰ, ρₒ, cₒ, liquidus::LinearLiquidus)

Solve the three-equation system for interface temperature, salinity, and melt rate.

The three equations are:
1. Heat balance: ``ρₒ cₒ αₕ u★ (Tₒ - Tᵦ) = ℰ q``
2. Salt balance: ``αₛ u★ (Sₒ - Sᵦ) = q (Sᵦ - Sᵢ)``
3. Freezing point: ``Tᵦ = Tₘ(Sᵦ)``

Returns `(Tᵦ, Sᵦ, q)` where q is the melt rate (positive for melting).
"""
@inline function solve_interface_conditions(Tₒ, Sₒ, Sᵢ, αₕ, αₛ, u★, ℰ, ρₒ, cₒ, liquidus::LinearLiquidus)
    # Liquidus: Tₘ(S) = Tₘ₀ - m * S where m = liquidus.slope
    # So Tₘ(S) = λ₁ * S + λ₂ with λ₁ = -m, λ₂ = Tₘ₀
    λ₁ = -liquidus.slope
    λ₂ = liquidus.freshwater_melting_temperature

    # Combine equations to get quadratic in Sᵦ:
    # From heat balance: q = ρₒ * cₒ * αₕ * u★ * (Tₒ - Tᵦ) / L
    # Substituting Tᵦ = λ₁ * Sᵦ + λ₂:
    #   q = ρₒ * cₒ * αₕ * u★ * (Tₒ - λ₁ * Sᵦ - λ₂) / L
    # Substituting into salt balance:
    #   αₛ * u★ * (Sₒ - Sᵦ) = q * (Sᵦ - Sᵢ)
    # This gives a quadratic: a * Sᵦ² + b * Sᵦ + c = 0

    # Coefficient: η = ρₒ * cₒ * αₕ * u★ / L
    η = ρₒ * cₒ * αₕ * u★ / ℰ
    αₛu★ = αₛ * u★

    # Quadratic coefficients
    a = λ₁ * η
    b = -αₛu★ - η * (Tₒ - λ₂ + λ₁ * Sᵢ)
    c = αₛu★ * Sₒ + η * (Tₒ - λ₂) * Sᵢ

    # Pre-compute reciprocal with zero check (MITgcm approach).
    # When a = 0 (e.g., u★ = 0), this avoids division by zero.
    ξ = ifelse(a == zero(a), zero(a), one(a) / (2a))

    # Solve quadratic: Sᵦ = (-b ± √Δ) / (2a)
    Δ = b^2 - 4a * c
    Δ = max(Δ, zero(Δ))

    # Try the root with -√Δ first (typically the physically meaningful one for a < 0)
    Sᵦ = (-b - sqrt(Δ)) * ξ

    # If this root yields negative salinity, use the other root (MITgcm approach)
    Sᵦ = ifelse(Sᵦ < zero(Sᵦ), (-b + sqrt(Δ)) * ξ, Sᵦ)

    # Interface temperature from liquidus
    Tᵦ = λ₁ * Sᵦ + λ₂

    # Melt rate from heat balance
    q = η * (Tₒ - Tᵦ)

    return Tᵦ, Sᵦ, q
end

#####
##### Show methods
#####

Base.summary(::IceBathHeatFlux{FT}) where FT = "IceBathHeatFlux{$FT}"
Base.summary(::ThreeEquationHeatFlux{FT}) where FT = "ThreeEquationHeatFlux{$FT}"

function Base.show(io::IO, flux::IceBathHeatFlux)
    print(io, summary(flux), '\n')
    print(io, "├── heat_transfer_coefficient: ", flux.heat_transfer_coefficient, '\n')
    print(io, "└── friction_velocity: ", summary(flux.friction_velocity))
end

function Base.show(io::IO, flux::ThreeEquationHeatFlux)
    print(io, summary(flux), '\n')
    print(io, "├── heat_transfer_coefficient: ", flux.heat_transfer_coefficient, '\n')
    print(io, "├── salt_transfer_coefficient: ", flux.salt_transfer_coefficient, '\n')
    print(io, "└── friction_velocity: ", summary(flux.friction_velocity))
end
