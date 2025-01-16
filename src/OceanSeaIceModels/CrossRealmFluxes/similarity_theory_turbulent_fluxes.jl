using Oceananigans.Utils: prettysummary
using Oceananigans.Grids: AbstractGrid

using Adapt
using Thermodynamics: Liquid
using SurfaceFluxes.Parameters: SurfaceFluxesParameters
using SurfaceFluxes.UniversalFunctions: BusingerParams, BusingerType

using Printf
using Thermodynamics: PhasePartition
using KernelAbstractions.Extras.LoopInfo: @unroll

using ..PrescribedAtmospheres: PrescribedAtmosphereThermodynamicsParameters

using Statistics: norm

import Thermodynamics as AtmosphericThermodynamics
import Thermodynamics.Parameters: molmass_ratio

#####
##### Bulk turbulent fluxes based on similarity theory
#####

struct SimilarityTheoryFluxes{FT, UF, R, B, T, V}
    gravitational_acceleration :: FT # parameter
    von_karman_constant :: FT        # parameter
    turbulent_prandtl_number :: FT   # parameter
    gustiness_parameter :: FT        # bulk velocity parameter
    stability_functions :: UF        # functions for turbulent fluxes
    roughness_lengths :: R           # parameterization for turbulent fluxes
    similarity_profile_type :: B     # similarity profile relating atmosphere to interface state
    interface_temperature_type :: T  # interface temperature either diagnostic or prescribed
    bulk_velocity :: V               # bulk velocity scale for turbulent fluxes
    tolerance :: FT                  # solver option
    maxiter :: Int                   # solver option
end

Adapt.adapt_structure(to, fluxes::SimilarityTheoryFluxes) = 
    SimilarityTheoryFluxes(adapt(to, fluxes.von_karman_constant),
                           adapt(to, fluxes.turbulent_prandtl_number),
                           adapt(to, fluxes.gustiness_parameter),
                           adapt(to, fluxes.stability_functions),
                           adapt(to, fluxes.roughness_lengths),
                           adapt(to, fluxes.similarity_profile_type),
                           adapt(to, fluxes.interface_temperature_type),
                           adapt(to, fluxes.bulk_velocity),
                           fluxes.tolerance,
                           fluxes.maxiter)

Base.summary(::SimilarityTheoryFluxes{FT}) where FT = "SimilarityTheoryFluxes{$FT}"

function Base.show(io::IO, fluxes::SimilarityTheoryFluxes)
    print(io, summary(fluxes), '\n',
          "├── von_karman_constant: ",        prettysummary(fluxes.von_karman_constant), '\n',
          "├── turbulent_prandtl_number: ",   prettysummary(fluxes.turbulent_prandtl_number), '\n',
          "├── gustiness_parameter: ",        prettysummary(fluxes.gustiness_parameter), '\n',
          "├── stability_functions: ",        summary(fluxes.stability_functions), '\n',
          "├── water_mole_fraction: ",        summary(fluxes.water_mole_fraction), '\n',
          "├── water_vapor_saturation: ",     summary(fluxes.water_vapor_saturation), '\n',
          "├── roughness_lengths: ",          summary(fluxes.roughness_lengths), '\n',
          "├── similarity_profile_type: ",    summary(fluxes.similarity_profile_type), '\n',
          "├── interface_temperature: ",      summary(fluxes.interface_temperature_type), '\n',
          "└── thermodynamics_parameters: ",  summary(fluxes.thermodynamics_parameters))
end

"""
    SimilarityTheoryFluxes(FT::DataType = Float64;
                           gravitational_acceleration = convert(FT, 9.81),
                           von_karman_constant = convert(FT, 0.4),
                           turbulent_prandtl_number = convert(FT, 1),
                           gustiness_parameter = convert(FT, 6.5),
                           stability_functions = default_stability_functions(FT),
                           roughness_lengths = default_roughness_lengths(FT),
                           similarity_profile_type = LogarithmicSimilarityProfile(),
                           interface_temperature_type = BulkTemperature(),
                           bulk_velocity = RelativeVelocity(),
                           tolerance = 1e-8,
                           maxiter = 100)

`SimilarityTheoryFluxes` contains parameters and settings to calculate
air-interface turbulent fluxes using Monin-Obukhov similarity theory.

Keyword Arguments
==================

- `von_karman_constant`: The von Karman constant. Default: 0.4.
- `turbulent_prandtl_number`: The turbulent Prandtl number. Default: 1.
- `gustiness_parameter`: The gustiness parameter that accounts for low wind speed areas. Default: 6.5.
- `stability_functions`: The stability functions. Default: `default_stability_functions(FT)` that follow the 
                         formulation of Edson et al. (2013).
- `roughness_lengths`: The roughness lengths used to calculate the characteristic scales for momentum, temperature and 
                       water vapor. Default: `default_roughness_lengths(FT)`, formulation taken from Edson et al (2013).
- `similarity_profile_type`: The type of similarity profile used to relate the atmospheric state to the 
                             interface fluxes / characteristic scales.
- `bulk_velocity`: The velocity used to calculate the characteristic scales. Default: `RelativeVelocity()` (difference between
                   atmospheric and interfaceic speed).
- `tolerance`: The tolerance for convergence. Default: 1e-8.
- `maxiter`: The maximum number of iterations. Default: 100.
"""
function SimilarityTheoryFluxes(FT::DataType = Float64;
                                gravitational_acceleration = 9.81,
                                von_karman_constant = 0.4,
                                turbulent_prandtl_number = 1,
                                gustiness_parameter = 6.5,
                                stability_functions = edson_stability_functions(FT),
                                roughness_lengths = default_roughness_lengths(FT),
                                similarity_profile_type = LogarithmicSimilarityProfile(),
                                interface_temperature_type = BulkTemperature(),
                                bulk_velocity = RelativeVelocity(),
                                tolerance = 1e-8,
                                maxiter = 100)

    return SimilarityTheoryFluxes(convert(FT, gravitational_acceleration),
                                  convert(FT, von_karman_constant),
                                  convert(FT, turbulent_prandtl_number),
                                  convert(FT, gustiness_parameter),
                                  stability_functions,
                                  roughness_lengths,
                                  similarity_profile_type,
                                  interface_temperature_type,
                                  bulk_velocity,
                                  convert(FT, tolerance), 
                                  maxiter)
end

#####
##### Similarity profile types
#####

"""
    LogarithmicSimilarityProfile()

Represent the classic Monin-Obukhov similarity profile, which finds that 

```math
ϕ(z) = Π(z) ϕ★ / ϰ
```

where ``ϰ`` is the Von Karman constant, ``ϕ★`` is the characteristic scale for ``ϕ``,
and ``Π`` is the "similarity profile",

```math
Π(h) = log(h / ℓ) - ψ(h / L) + ψ(ℓ / L)
```

which is a logarithmic profile adjusted by the stability function ``ψ`` and dependent on
the Monin-Obukhov length ``L`` and the roughness length ``ℓ``.
"""
struct LogarithmicSimilarityProfile end
struct COARELogarithmicSimilarityProfile end

@inline similarity_profile(::LogarithmicSimilarityProfile, ψ, h, ℓ, L) =
    log(h / ℓ) - ψ(h / L) + ψ(ℓ / L)

@inline similarity_profile(::COARELogarithmicSimilarityProfile, ψ, h, ℓ, L) =
    log(h / ℓ) - ψ(h / L)

# Iterating condition for the characteristic scales solvers
@inline function _iterating(Ψⁿ, Ψ⁻, iteration, maxiter, tolerance)
    hasnt_started = iteration == 0
    reached_maxiter = iteration ≥ maxiter
    drift = abs(Ψⁿ.u★ - Ψ⁻.u★) + abs(Ψⁿ.θ★ - Ψ⁻.θ★) + abs(Ψⁿ.q★ - Ψ⁻.q★)
    converged = drift < tolerance
    return !(converged | reached_maxiter) | hasnt_started
end

@inline function compute_interface_state(turbulent_flux_formulation,
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
    maxiter = turbulent_flux_formulation.maxiter
    tolerance = turbulent_flux_formulation.tolerance

    while _iterating(Ψₛⁿ, Ψₛ⁻, iteration, maxiter, tolerance)
        Ψₛ⁻ = Ψₛⁿ
        Ψₛⁿ = iterate_interface_state(turbulent_flux_formulation,
                                      Ψₛ⁻, Ψₐ, Ψᵢ,
                                      downwelling_radiation,
                                      interface_properties,
                                      atmosphere_properties,
                                      interior_properties)
        iteration += 1
    end

    return Ψₛⁿ

end

@inline function iterate_interface_state(turbulent_flux_formulation,
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
    g = turbulent_flux_formulation.gravitational_acceleration
    cₚ = interior_properties.heat_capacity
    Δθ = Tₐ - Tₛ + g / cₚ * Δh

    # Recompute interface thermodynamic state with new temperature and specific humidity
    𝒬ₛ = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, 𝒬ₐ.p, Tₛ, qₛ)

    # "initial" scales because we will recompute them
    u★ = approximate_interface_state.u★
    θ★ = approximate_interface_state.θ★
    q★ = approximate_interface_state.q★

    # Similarity functions from Edson et al. (2013)
    ψu = turbulent_flux_formulation.stability_functions.momentum
    ψθ = turbulent_flux_formulation.stability_functions.temperature
    ψq = turbulent_flux_formulation.stability_functions.water_vapor

    # Extract roughness lengths
    ℓu = turbulent_flux_formulation.roughness_lengths.momentum
    ℓθ = turbulent_flux_formulation.roughness_lengths.temperature
    ℓq = turbulent_flux_formulation.roughness_lengths.water_vapor
    β = turbulent_flux_formulation.gustiness_parameter

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    b★ = buoyancy_scale(θ★, q★, 𝒬ₛ, ℂₐ, g)

    # Monin-Obhukov characteristic length scale and non-dimensional height
    ϰ = turbulent_flux_formulation.von_karman_constant
    L★ = ifelse(b★ == 0, zero(b★), - u★^2 / (ϰ * b★))

    # Compute roughness length scales
    ℓu₀ = roughness_length(ℓu, u★, 𝒬ₛ, ℂₐ)
    ℓq₀ = roughness_length(ℓq, ℓu₀, u★, 𝒬ₛ, ℂₐ)
    ℓθ₀ = roughness_length(ℓθ, ℓu₀, u★, 𝒬ₛ, ℂₐ)

    # Transfer coefficients at height `h`
    profile_type = turbulent_flux_formulation.similarity_profile_type
    χu = ϰ / similarity_profile(profile_type, ψu, Δh, ℓu₀, L★)
    χθ = ϰ / similarity_profile(profile_type, ψθ, Δh, ℓθ₀, L★)
    χq = ϰ / similarity_profile(profile_type, ψq, Δh, ℓq₀, L★)

    # Buoyancy flux characteristic scale for gustiness (Edson 2013)
    h_bℓ = atmosphere_state.h_bℓ
    Jᵇ = - u★ * b★
    Uᴳ = β * cbrt(Jᵇ * h_bℓ)

    # New velocity difference accounting for gustiness
    Δu, Δv = velocity_difference(turbulent_flux_formulation.bulk_velocity, atmosphere_state, approximate_interface_state)
    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳ^2)

    #=
    Pr = turbulent_flux_formulation.turbulent_prandtl_number
    χθ = χθ / Pr
    χq = χq / Pr
    =#

    # Recompute 
    u★ = χu * ΔU
    θ★ = χθ * Δθ
    q★ = χq * Δq

    u = approximate_interface_state.u
    v = approximate_interface_state.v
    S = approximate_interface_state.S

    return InterfaceState(u★, θ★, q★, u, v, Tₛ, S, convert(FT, qₛ))
end

""" The exchange fluxes depend on the atmosphere velocity but not the interface velocity """
struct WindVelocity end

""" The exchange fluxes depend on the relative velocity between the atmosphere and the interface """
struct RelativeVelocity end

"""
    buoyancy_scale(θ★, q★, 𝒬, ℂ, g)

Return the characteristic buoyancy scale `b★` associated with
the characteristic temperature `θ★`, specific humidity scale `q★`,
near-interface atmospheric thermodynamic state `𝒬`, thermodynamic
parameters `ℂ`, and gravitational acceleration `g`.

The buoyancy scale is defined in terms of the interface buoyancy flux,

```math
u★ b★ ≡ w′b′,
```

where `u★` is the friction velocity.
Using the definition of buoyancy for non-condensing air, we find that

```math
b★ = g / 𝒯ₐ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₐ * q★),
```
where ``𝒯ₐ`` is the virtual temperature of the atmosphere near the interface,
and ``δ = Rᵥ / R_d - 1``, where ``Rᵥ`` is the molar mass of water vapor and
``R_d`` is the molar mass of dry air.

Note that the Monin-Obukhov characteristic length scale is defined
in terms of `b★` and additionally the Von Karman constant `ϰ`,

```math
L★ = - u★² / ϰ b★ .
```
"""
@inline function buoyancy_scale(θ★, q★, 𝒬, ℂ, g)
    𝒯ₐ = AtmosphericThermodynamics.virtual_temperature(ℂ, 𝒬)
    qₐ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬)
    ε  = AtmosphericThermodynamics.Parameters.molmass_ratio(ℂ)
    δ  = ε - 1 # typically equal to 0.608

    b★ = g / 𝒯ₐ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₐ * q★)

    return b★
end

@inline function velocity_difference(::RelativeVelocity, 𝒰₁, 𝒰₀)
    Δu = 𝒰₁.u - 𝒰₀.u
    Δv = 𝒰₁.v - 𝒰₀.v
    return Δu, Δv
end

@inline velocity_difference(::WindVelocity, 𝒰₁, 𝒰₀) = 𝒰₁.u, 𝒰₁.v

