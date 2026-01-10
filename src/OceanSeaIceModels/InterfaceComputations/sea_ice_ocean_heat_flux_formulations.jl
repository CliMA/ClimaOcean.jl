using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

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
1. Heat balance: ``Q_o = \\rho_o c_o \\alpha_h u_* (T - T_b)``
2. Salt balance: ``(S_b - S_i) G_h = \\alpha_s u_* (S - S_b)``
3. Freezing point: ``T_b = T_m(S_b)``

where ``T_b`` and ``S_b`` are the interface temperature and salinity,
``\\alpha_h`` and ``\\alpha_s`` are heat and salt transfer coefficients,
and ``G_h`` is the ice thermodynamic tendency (growth/melt rate).

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
    compute_interface_heat_flux(flux, i, j, Tᵢ, Sᵢ, Tₒ, Sₒ, Sⁱ, ℵ, Gₕ, Nz, liquidus, ρₒ, cₒ, τx, τy)

Compute the heat flux at the sea ice-ocean interface for grid point `(i, j)`.
Returns the heat flux ``Q`` where ``Q > 0`` means heat flux from ocean to ice (ocean cooling).
"""
@inline function compute_interface_heat_flux(flux::IceBathHeatFlux, i, j,
                                              Tᵢ, Sᵢ, Tₒ, Sₒ, Sⁱ, ℵ, Gₕ, Nz,
                                              liquidus, ρₒ, cₒ, τx, τy)
    @inbounds begin
        Tᴺ  = Tᵢ[i, j, 1]
        Sᴺ  = Sᵢ[i, j, 1]
        ℵᵢⱼ = ℵ[i, j, 1]
    end

    # Interface temperature is at the freezing point
    Tₘ = melting_temperature(liquidus, Sᴺ)

    αₕ = flux.heat_transfer_coefficient
    u★ = get_friction_velocity(flux.friction_velocity, i, j, τx, τy, ρₒ)

    # Q > 0 means heat flux from ocean to ice (ocean cooling)
    return ρₒ * cₒ * αₕ * u★ * (Tᴺ - Tₘ) * ℵᵢⱼ
end

@inline function compute_interface_heat_flux(flux::ThreeEquationHeatFlux, i, j,
                                              Tᵢ, Sᵢ, Tₒ, Sₒ, Sⁱ, ℵ, Gₕ, Nz,
                                              liquidus, ρₒ, cₒ, τx, τy)
    @inbounds begin
        Tᴺ   = Tₒ[i, j, Nz]   # Ocean surface temperature
        Sᴺ   = Sₒ[i, j, Nz]   # Ocean surface salinity
        Sᵢᶜᵉ = Sⁱ[i, j, 1]    # Ice salinity
        ℵᵢⱼ  = ℵ[i, j, 1]     # Ice concentration
        gₕ   = Gₕ[i, j, 1]     # Ice growth rate
    end

    αₕ = flux.heat_transfer_coefficient
    αₛ = flux.salt_transfer_coefficient
    u★ = get_friction_velocity(flux.friction_velocity, i, j, τx, τy, ρₒ)

    # Solve for interface temperature and salinity
    Tᵦ, Sᵦ = solve_interface_conditions(Tᴺ, Sᴺ, Sᵢᶜᵉ, gₕ, αₕ, αₛ, u★, liquidus)

    # Store interface values
    @inbounds Tᵢ[i, j, 1] = Tᵦ
    @inbounds Sᵢ[i, j, 1] = Sᵦ

    # Q > 0 means heat flux from ocean to ice (ocean cooling)
    return ρₒ * cₒ * αₕ * u★ * (Tᴺ - Tᵦ) * ℵᵢⱼ
end

"""
    solve_interface_conditions(Tₒ, Sₒ, Sᵢ, Gₕ, αₕ, αₛ, u★, liquidus)

Solve the three-equation system for interface temperature and salinity.

The system consists of:
1. Heat balance at interface
2. Salt balance at interface: ``(S_b - S_i) G_h = \\alpha_s u_* (S_o - S_b)``
3. Freezing point constraint: ``T_b = T_m(S_b)``

For a linear liquidus ``T_m(S) = \\lambda_1 S + \\lambda_2``, an analytical solution exists.
"""
@inline function solve_interface_conditions(Tₒ, Sₒ, Sᵢ, Gₕ, αₕ, αₛ, u★, liquidus)
    # For a linear liquidus: Tₘ(S) = λ₁ * S + λ₂
    # The salt balance gives: (Sᵦ - Sᵢ) * Gₕ = αₛ * u★ * (Sₒ - Sᵦ)
    # Solving for Sᵦ:
    #   Sᵦ * Gₕ - Sᵢ * Gₕ = αₛ * u★ * Sₒ - αₛ * u★ * Sᵦ
    #   Sᵦ * (Gₕ + αₛ * u★) = αₛ * u★ * Sₒ + Sᵢ * Gₕ
    #   Sᵦ = (αₛ * u★ * Sₒ + Sᵢ * Gₕ) / (Gₕ + αₛ * u★)

    αₛu★ = αₛ * u★

    # Handle the case when Gₕ + αₛu★ ≈ 0 (no ice growth and no turbulent transfer)
    denominator = Gₕ + αₛu★

    # If denominator is very small, fall back to ocean salinity
    Sᵦ = ifelse(abs(denominator) < eps(typeof(denominator)), Sₒ,
                (αₛu★ * Sₒ + Sᵢ * Gₕ) / denominator)

    # Interface temperature from liquidus
    Tᵦ = melting_temperature(liquidus, Sᵦ)

    return Tᵦ, Sᵦ
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
