using Oceananigans.Utils: prettysummary
using Oceananigans.Grids: AbstractGrid
using Oceananigans.BuoyancyFormulations: g_Earth

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
##### These are more general properties
#####

""" The exchange fluxes depend on the atmosphere velocity but not the interface velocity """
struct WindVelocity end

""" The exchange fluxes depend on the relative velocity between the atmosphere and the interface """
struct RelativeVelocity end


#####
##### Bulk turbulent fluxes based on similarity theory
#####

struct SimilarityTheoryFluxes{FT, UF, R, B, V}
    gravitational_acceleration :: FT # parameter
    von_karman_constant :: FT        # parameter
    turbulent_prandtl_number :: FT   # parameter
    gustiness_parameter :: FT        # bulk velocity parameter
    stability_functions :: UF        # functions for turbulent fluxes
    roughness_lengths :: R           # parameterization for turbulent fluxes
    similarity_form :: B             # similarity profile relating atmosphere to interface state
    bulk_velocity :: V               # bulk velocity scale for turbulent fluxes
    solver_tolerance :: FT           # solver option
    solver_maxiter :: Int            # solver option
end

Adapt.adapt_structure(to, fluxes::SimilarityTheoryFluxes) = 
    SimilarityTheoryFluxes(adapt(to, fluxes.gravitational_acceleration),
                           adapt(to, fluxes.von_karman_constant),
                           adapt(to, fluxes.turbulent_prandtl_number),
                           adapt(to, fluxes.gustiness_parameter),
                           adapt(to, fluxes.stability_functions),
                           adapt(to, fluxes.roughness_lengths),
                           adapt(to, fluxes.similarity_form),
                           adapt(to, fluxes.bulk_velocity),
                           fluxes.solver_tolerance,
                           fluxes.solver_maxiter)

Base.summary(::SimilarityTheoryFluxes{FT}) where FT = "SimilarityTheoryFluxes{$FT}"

function Base.show(io::IO, fluxes::SimilarityTheoryFluxes)
    print(io, summary(fluxes), '\n',
          "├── gravitational_acceleration: ", prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",        prettysummary(fluxes.von_karman_constant), '\n',
          "├── turbulent_prandtl_number: ",   prettysummary(fluxes.turbulent_prandtl_number), '\n',
          "├── gustiness_parameter: ",        prettysummary(fluxes.gustiness_parameter), '\n',
          "├── stability_functions: ",        summary(fluxes.stability_functions), '\n',
          "├── roughness_lengths: ",          summary(fluxes.roughness_lengths), '\n',
          "├── bulk_velocity: ",              summary(fluxes.bulk_velocity), '\n',
          "├── similarity_form: ",            summary(fluxes.similarity_form), '\n',
          "├── solver_tolerance: ",           summary(fluxes.solver_tolerance), '\n',
          "└── solver_maxiter: ",             summary(fluxes.solver_maxiter))
end

"""
    SimilarityTheoryFluxes(FT::DataType = Float64;
                           gravitational_acceleration = convert(FT, 9.81),
                           von_karman_constant = convert(FT, 0.4),
                           turbulent_prandtl_number = convert(FT, 1),
                           gustiness_parameter = convert(FT, 6.5),
                           stability_functions = default_stability_functions(FT),
                           roughness_lengths = default_roughness_lengths(FT),
                           similarity_form = LogarithmicSimilarityProfile(),
                           bulk_velocity = RelativeVelocity(),
                           solver_tolerance = 1e-8,
                           solver_maxiter = 100)

`SimilarityTheoryFluxes` contains parameters and settings to calculate
air-interface turbulent fluxes using Monin-Obukhov similarity theory.

Keyword Arguments
==================

- `gravitational_acceleration`: Gravitational acceleration.
- `von_karman_constant`: The von Karman constant. Default: 0.4.
- `turbulent_prandtl_number`: The turbulent Prandtl number. Default: 1.
- `gustiness_parameter`: The gustiness parameter that accounts for low wind speed areas. Default: 6.5.
- `stability_functions`: The stability functions. Default: `default_stability_functions(FT)` that follow the 
                         formulation of Edson et al. (2013).
- `roughness_lengths`: The roughness lengths used to calculate the characteristic scales for momentum, temperature and 
                       water vapor. Default: `default_roughness_lengths(FT)`, formulation taken from Edson et al (2013).
- `similarity_form`: The type of similarity profile used to relate the atmospheric state to the 
                             interface fluxes / characteristic scales.
- `bulk_velocity`: The velocity used to calculate the characteristic scales. Default: `RelativeVelocity()` (difference between
                   atmospheric and interface speed).
- `solver_tolerance`: The tolerance for convergence. Default: 1e-8.
- `solver_maxiter`: The maximum number of iterations. Default: 100.
"""
function SimilarityTheoryFluxes(FT::DataType = Oceananigans.defaults.FloatType;
                                gravitational_acceleration = g_Earth,
                                von_karman_constant = 0.4,
                                turbulent_prandtl_number = 1,
                                gustiness_parameter = 6.5,
                                stability_functions = edson_stability_functions(FT),
                                roughness_lengths = default_roughness_lengths(FT),
                                similarity_form = LogarithmicSimilarityProfile(),
                                bulk_velocity = RelativeVelocity(),
                                solver_tolerance = 1e-8,
                                solver_maxiter = 100)

    return SimilarityTheoryFluxes(convert(FT, gravitational_acceleration),
                                  convert(FT, von_karman_constant),
                                  convert(FT, turbulent_prandtl_number),
                                  convert(FT, gustiness_parameter),
                                  stability_functions,
                                  roughness_lengths,
                                  similarity_form,
                                  bulk_velocity,
                                  convert(FT, solver_tolerance), 
                                  solver_maxiter)
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

function iterate_interface_fluxes(flux_formulation::SimilarityTheoryFluxes,
                                  Tₛ, qₛ, Δθ, Δq, Δh,
                                  approximate_interface_state,
                                  atmosphere_state,
                                  atmosphere_properties)

    ℂₐ = atmosphere_properties.thermodynamics_parameters
    𝒬ₐ = atmosphere_state.𝒬

    # "initial" scales because we will recompute them
    u★ = approximate_interface_state.u★
    θ★ = approximate_interface_state.θ★
    q★ = approximate_interface_state.q★

    # Similarity functions from Edson et al. (2013)
    ψu = flux_formulation.stability_functions.momentum
    ψθ = flux_formulation.stability_functions.temperature
    ψq = flux_formulation.stability_functions.water_vapor

    # Extract roughness lengths
    ℓu = flux_formulation.roughness_lengths.momentum
    ℓθ = flux_formulation.roughness_lengths.temperature
    ℓq = flux_formulation.roughness_lengths.water_vapor
    β = flux_formulation.gustiness_parameter

    # Compute surface thermodynamic state
    𝒬ₛ = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, 𝒬ₐ.p, Tₛ, qₛ)

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    g = flux_formulation.gravitational_acceleration
    b★ = buoyancy_scale(θ★, q★, 𝒬ₛ, ℂₐ, g)

    # Monin-Obhukov characteristic length scale and non-dimensional height
    ϰ = flux_formulation.von_karman_constant
    L★ = ifelse(b★ == 0, Inf, - u★^2 / (ϰ * b★))

    # Compute roughness length scales
    ℓu₀ = roughness_length(ℓu, u★, 𝒬ₛ, ℂₐ)
    ℓq₀ = roughness_length(ℓq, ℓu₀, u★, 𝒬ₛ, ℂₐ)
    ℓθ₀ = roughness_length(ℓθ, ℓu₀, u★, 𝒬ₛ, ℂₐ)

    # Transfer coefficients at height `h`
    form = flux_formulation.similarity_form
    χu = ϰ / similarity_profile(form, ψu, Δh, ℓu₀, L★)
    χθ = ϰ / similarity_profile(form, ψθ, Δh, ℓθ₀, L★)
    χq = ϰ / similarity_profile(form, ψq, Δh, ℓq₀, L★)

    #=
    Pr = flux_formulation.turbulent_prandtl_number
    χθ = χθ / Pr
    χq = χq / Pr
    =#
    
    # Buoyancy flux characteristic scale for gustiness (Edson 2013)
    h_bℓ = atmosphere_state.h_bℓ
    Jᵇ = - u★ * b★
    Uᴳ = β * cbrt(Jᵇ * h_bℓ)

    # New velocity difference accounting for gustiness
    Δu, Δv = velocity_difference(flux_formulation.bulk_velocity,
                                 atmosphere_state,
                                 approximate_interface_state)

    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳ^2)

    # Recompute 
    u★ = χu * ΔU
    θ★ = χθ * Δθ
    q★ = χq * Δq

    return u★, θ★, q★
end

"""
    buoyancy_scale(θ★, q★, 𝒬, ℂ, g)

Return the characteristic buoyancy scale `b★` associated with
the characteristic temperature `θ★`, specific humidity scale `q★`,
surface thermodynamic state `𝒬`, thermodynamic
parameters `ℂ`, and gravitational acceleration `g`.

The buoyancy scale is defined in terms of the interface buoyancy flux,

```math
u★ b★ ≡ w′b′,
```

where `u★` is the friction velocity.
Using the definition of buoyancy for non-condensing air, we find that

```math
b★ = g / 𝒯ₛ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₛ * q★),
```
where ``𝒯ₐ`` is the virtual temperature at the surface,
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

import Statistics

#####
##### Struct that represents a 3-tuple of momentum, heat, and water vapor
#####

struct SimilarityScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

"""
    EdsonMomentumStabilityFunction{FT}

A struct representing the momentum stability function detailed in Edson et al (2013).
The formulation hinges on the definition of three different functions:
one for stable atmospheric conditions ``(ζ > 0)``, named ``ψₛ`` and two for unstable conditions,
named ``ψᵤ₁`` and ``ψᵤ₂``.
These stability functions are obtained by regression to experimental data.

The stability parameter for stable atmospheric conditions is defined as
```math
dζ = min(ζmax, Aˢζ)
ψₛ = - (Bˢ ζ + Cˢ ( ζ - Dˢ ) ) exp( - dζ) - Cˢ Dˢ 
```

While the stability parameter for unstable atmospheric conditions is calculated
as a function of the two individual stability functions as follows

```math
fᵤ₁ = √√(1 - Aᵘζ)
ψᵤ₁ = Bᵘ / 2 ⋅ log((1 + fᵤ₁ + fᵤ₁² + fᵤ₁³) / Bᵘ) - √Bᵘ atan(fᵤ₁) - Cᵘ

fᵤ₂ = ∛(1 - Dᵘζ)
ψᵤ₂ = Eᵘ / 2 ⋅ log((1 + fᵤ₂ + fᵤ₂²) / Eᵘ) - √Eᵘ atan( (1 + 2fᵤ₂) / √Eᵘ) + Fᵘ

f  = ζ² / (1 + ζ²)
ψᵤ = (1 - f) ψᵤ₁ + f ψᵤ₂  
```

The superscripts ``ˢ`` and ``ᵘ`` indicate if the parameter applies to the 
stability function for _stable_ or _unstable_ atmospheric conditions, respectively.
"""
@kwdef struct EdsonMomentumStabilityFunction{FT}
    ζmax :: FT = 50.0
    Aˢ   :: FT = 0.35
    Bˢ   :: FT = 0.7
    Cˢ   :: FT = 0.75
    Dˢ   :: FT = 5/0.35
    Aᵘ   :: FT = 15.0
    Bᵘ   :: FT = 2.0
    Cᵘ   :: FT = π/2
    Dᵘ   :: FT = 10.15
    Eᵘ   :: FT = 3.0
    Fᵘ   :: FT = π / sqrt(3)
end

@inline function (ψ::EdsonMomentumStabilityFunction)(ζ)
    ζmax = ψ.ζmax
    Aˢ   = ψ.Aˢ  
    Bˢ   = ψ.Bˢ  
    Cˢ   = ψ.Cˢ  
    Dˢ   = ψ.Dˢ  
    Aᵘ   = ψ.Aᵘ  
    Bᵘ   = ψ.Bᵘ  
    Cᵘ   = ψ.Cᵘ  
    Dᵘ   = ψ.Dᵘ  
    Eᵘ   = ψ.Eᵘ  
    Fᵘ   = ψ.Fᵘ  

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(ζmax, Aˢ * ζ⁺)

    # Stability parameter for _stable_ atmospheric conditions
    ψₛ = - (Bˢ * ζ⁺ + Cˢ * (ζ⁺ - Dˢ)) * exp(- dζ) - Cˢ * Dˢ
        
    # Stability parameter for _unstable_ atmospheric conditions
    fᵤ₁ = sqrt(sqrt(1 - Aᵘ * ζ⁻))
    ψᵤ₁ = Bᵘ * log((1 + fᵤ₁) / Bᵘ) + log((1 + fᵤ₁^2) / Bᵘ) - Bᵘ * atan(fᵤ₁) + Cᵘ
        
    fᵤ₂ = cbrt(1 - Dᵘ * ζ⁻)
    ψᵤ₂ = Eᵘ / 2 * log((1 + fᵤ₂ + fᵤ₂^2) / Eᵘ) - sqrt(Eᵘ) * atan( (1 + 2fᵤ₂) / sqrt(Eᵘ)) + Fᵘ
        
    f  = ζ⁻^2 / (1 + ζ⁻^2)
    ψᵤ = (1 - f) * ψᵤ₁ + f * ψᵤ₂  
        
    return ifelse(ζ < 0, ψᵤ, ψₛ)
end

"""
    EdsonScalarStabilityFunction{FT}

A struct representing the scalar stability function detailed in Edson et al (2013).
The formulation hinges on the definition of three different functions:
one for stable atmospheric conditions ``(ζ > 0)``, named ``ψₛ`` and two for unstable conditions,
named ``ψᵤ₁`` and ``ψᵤ₂``.
These stability functions are obtained by regression to experimental data.

The stability parameter for stable atmospheric conditions is defined as
```math
dζ = min(ζmax, Aˢζ)
ψₛ = - (1 + Bˢ ζ) ^ Cₛ - Bˢ ( ζ - Dˢ ) * exp( - dζ) - Eˢ
```

While the stability parameter for unstable atmospheric conditions is calculated
as a function of the two individual stability functions as follows
```math
fᵤ₁ = √(1 - Aᵘζ)
ψᵤ₁ = Bᵘ ⋅ log((1 + fᵤ₁) / Bᵘ) + Cᵤ

fᵤ₂ = ∛(1 - Dᵘζ)
ψᵤ₂ = Eᵘ / 2 ⋅ log((1 + fᵤ₂ + fᵤ₂²) / Eᵘ) - √Eᵘ atan( (1 + 2fᵤ₂) / √Eᵘ) + Fᵘ

f  = ζ² / (1 + ζ²)
ψᵤ = (1 - f) ψᵤ₁ + f ψᵤ₂  
```

The superscripts ``ˢ`` and ``ᵘ`` indicate if the parameter applies to the 
stability function for _stable_ or _unstable_ atmospheric conditions, respectively.
"""
@kwdef struct EdsonScalarStabilityFunction{FT}
    ζmax :: FT = 50.0
    Aˢ   :: FT = 0.35
    Bˢ   :: FT = 2/3
    Cˢ   :: FT = 3/2
    Dˢ   :: FT = 14.28
    Eˢ   :: FT = 8.525
    Aᵘ   :: FT = 15.0
    Bᵘ   :: FT = 2.0
    Cᵘ   :: FT = 0.0
    Dᵘ   :: FT = 34.15
    Eᵘ   :: FT = 3.0
    Fᵘ   :: FT = π / sqrt(3)
end

@inline function (ψ::EdsonScalarStabilityFunction)(ζ)
    ζmax = ψ.ζmax
    Aˢ   = ψ.Aˢ  
    Bˢ   = ψ.Bˢ  
    Cˢ   = ψ.Cˢ  
    Dˢ   = ψ.Dˢ  
    Eˢ   = ψ.Eˢ  
    Aᵘ   = ψ.Aᵘ  
    Bᵘ   = ψ.Bᵘ  
    Cᵘ   = ψ.Cᵘ  
    Dᵘ   = ψ.Dᵘ  
    Eᵘ   = ψ.Eᵘ  
    Fᵘ   = ψ.Fᵘ  

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(ζmax, Aˢ * ζ⁺)

    # stability function for stable atmospheric conditions 
    ψₛ = - (1 + Bˢ * ζ⁺)^Cˢ - Bˢ * (ζ⁺ - Dˢ) * exp(-dζ) - Eˢ

    # Stability parameter for _unstable_ atmospheric conditions
    fᵤ₁ = sqrt(1 - Aᵘ * ζ⁻)
    ψᵤ₁ = Bᵘ * log((1 + fᵤ₁) / Bᵘ) + Cᵘ

    fᵤ₂ = cbrt(1 - Dᵘ * ζ⁻)
    ψᵤ₂ = Eᵘ / 2 * log((1 + fᵤ₂ + fᵤ₂^2) / Eᵘ) - sqrt(Eᵘ) * atan((1 + 2fᵤ₂) / sqrt(Eᵘ)) + Fᵘ

    f  = ζ⁻^2 / (1 + ζ⁻^2)
    ψᵤ = (1 - f) * ψᵤ₁ + f * ψᵤ₂  

    return ifelse(ζ < 0, ψᵤ, ψₛ)
end

# Edson et al. (2013)
function edson_stability_functions(FT=Oceananigans.defaults.FloatType)
    ψu = EdsonMomentumStabilityFunction{FT}()
    ψc = EdsonScalarStabilityFunction{FT}()
    return SimilarityScales(ψu, ψc, ψc)
end

#####
##### From Grachev et al 2007, for stable boundary layers
#####

@kwdef struct ShebaMomentumStabilityFunction{FT}
    a :: FT = 6.5
    b :: FT = 1.3
end

# @inline (ψ::ShebaMomentumStabilityFunction)(ζ) = 1 + ψ.a * ζ * cbrt(1 + ζ) / (ψ.b + ζ)
@inline function (Ψ::ShebaMomentumStabilityFunction)(ζ)
    a = Ψ.a
    b = Ψ.b
    ζ⁺ = max(zero(ζ), ζ)
    z = cbrt(1 + ζ⁺)
    B = cbrt((1 - b) / b)

    rt3 = sqrt(3)
    Ψ₁ = - 3 * a * (z - 1) / b
    Ψ₂ = a * B / 2b * (2 * log((z + B) / (1 + B))
                       - log((z^2 - B * z + B^2) / (1 - B + B^2))
                       + 2 * rt3 * (atan((2z - B) / (rt3 * B)) - atan((2 - B) / (rt3 * B))))

    return Ψ₁ + Ψ₂
end

@kwdef struct ShebaScalarStabilityFunction{FT}
    a :: FT = 5.0
    b :: FT = 5.0
    c :: FT = 3.0
end

@inline function (Ψ::ShebaScalarStabilityFunction)(ζ)
    a = Ψ.a
    b = Ψ.b
    c = Ψ.c
    B = sqrt(c^2 - 4)
    ζ⁺ = max(zero(ζ), ζ)

    Ψ₁ = - b/2 * log(1 + c * ζ⁺ + ζ⁺^2)
    Ψ₂ = (b * c / 2B - a / B) * (log((2ζ⁺ + c - B) / (2ζ⁺ + c + B))
                                 + log((c - B) / (c + B)))

    return Ψ₁ + Ψ₂
end

#####
##### From Paulson 1970 for unstable boundary layers
####

@kwdef struct PaulsonMomentumStabilityFunction{FT}
    a :: FT = 16.0
    b :: FT = π/2
end

@inline function (Ψ::PaulsonMomentumStabilityFunction)(ζ)
    a = Ψ.a
    b = Ψ.b
    ζ⁻ = min(zero(ζ), ζ)
    z = sqrt(sqrt((1 - a * ζ⁻)))

    Ψ₁ = 2 * log((1 + z) / 2)
    Ψ₂ = log((1 + z^2) / 2)
    Ψ₃ = - 2 * atan(z)

    return Ψ₁ + Ψ₂ + Ψ₃ + b
end

@kwdef struct PaulsonScalarStabilityFunction{FT}
    a :: FT = 16.0
end

@inline function (Ψ::PaulsonScalarStabilityFunction)(ζ)
    a = Ψ.a
    ζ⁻ = min(zero(ζ), ζ)
    z = sqrt(sqrt((1 - a * ζ⁻)))
    return 2 * log((1 + z^2) / 2)
end

struct SplitStabilityFunction{S, U}
    stable :: S
    unstable :: U
end

@inline function (Ψ::SplitStabilityFunction)(ζ)
    Ψ_stable = Ψ.stable(ζ)
    Ψ_unstable = Ψ.unstable(ζ)
    stable = ζ > 0
    return ifelse(stable, Ψ_stable, Ψ_unstable)
end

function atmosphere_sea_ice_stability_functions(FT=Oceananigans.defaults.FloatType)
    stable_momentum = PaulsonMomentumStabilityFunction{FT}()
    unstable_momentum = ShebaMomentumStabilityFunction{FT}()
    momentum = SplitStabilityFunction(stable_momentum, unstable_momentum)

    stable_scalar = PaulsonScalarStabilityFunction{FT}()
    unstable_scalar = ShebaScalarStabilityFunction{FT}()
    scalar = SplitStabilityFunction(stable_scalar, unstable_scalar)

    return SimilarityScales(momentum, scalar, scalar)
end

