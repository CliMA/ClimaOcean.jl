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

import SurfaceFluxes.Parameters:
    thermodynamics_params,
    uf_params,
    von_karman_const,
    grav

#####
##### Bulk turbulent fluxes based on similarity theory
#####

struct SimilarityTheoryTurbulentFluxes{FT, UF, TP, S, W, R, B, T, V, F}
    gravitational_acceleration :: FT # parameter
    von_karman_constant :: FT        # parameter
    turbulent_prandtl_number :: FT   # parameter
    gustiness_parameter :: FT        # bulk velocity parameter
    stability_functions :: UF        # functions for turbulent fluxes
    thermodynamics_parameters :: TP  # parameter group
    water_vapor_saturation :: S      # model for computing the saturation water vapor mass
    water_mole_fraction :: W         # mole fraction of H₂O in seawater
    roughness_lengths :: R           # parameterization for turbulent fluxes
    similarity_profile_type :: B     # similarity profile relating atmosphere to surface state
    surface_temperature_type :: T    # surface temperature either diagnostic or prescribed
    bulk_velocity :: V               # bulk velocity scale for turbulent fluxes
    tolerance :: FT                  # solver option
    maxiter :: Int                   # solver option
    fields :: F                      # fields that store turbulent fluxes
end

const STTF = SimilarityTheoryTurbulentFluxes
@inline thermodynamics_params(fluxes::STTF) = fluxes.thermodynamics_parameters
@inline uf_params(fluxes::STTF)             = fluxes.stability_functions
@inline von_karman_const(fluxes::STTF)      = fluxes.von_karman_constant
@inline grav(fluxes::STTF)                  = fluxes.gravitational_acceleration
@inline molmass_ratio(fluxes::STTF)         = molmass_ratio(fluxes.thermodynamics_parameters)

Adapt.adapt_structure(to, fluxes::STTF) = SimilarityTheoryTurbulentFluxes(adapt(to, fluxes.gravitational_acceleration),
                                                                          adapt(to, fluxes.von_karman_constant),
                                                                          adapt(to, fluxes.turbulent_prandtl_number),
                                                                          adapt(to, fluxes.gustiness_parameter),
                                                                          adapt(to, fluxes.stability_functions),
                                                                          adapt(to, fluxes.thermodynamics_parameters),
                                                                          adapt(to, fluxes.water_vapor_saturation),
                                                                          adapt(to, fluxes.water_mole_fraction),
                                                                          adapt(to, fluxes.roughness_lengths),
                                                                          adapt(to, fluxes.similarity_profile_type),
                                                                          adapt(to, fluxes.surface_temperature_type),
                                                                          adapt(to, fluxes.bulk_velocity),
                                                                          fluxes.tolerance,
                                                                          fluxes.maxiter,
                                                                          adapt(to, fluxes.fields))

Base.summary(::SimilarityTheoryTurbulentFluxes{FT}) where FT = "SimilarityTheoryTurbulentFluxes{$FT}"

struct ClasiusClapyeronSaturation end
 
@inline function water_saturation_specific_humidity(::ClasiusClapyeronSaturation, ℂₐ, ρₛ, Tₛ)
    FT = eltype(ℂₐ)
    p★ = AtmosphericThermodynamics.saturation_vapor_pressure(ℂₐ, convert(FT, Tₛ), Liquid())
    q★ = AtmosphericThermodynamics.q_vap_saturation_from_density(ℂₐ, convert(FT, Tₛ), ρₛ, p★)
    return q★
end

function Base.show(io::IO, fluxes::SimilarityTheoryTurbulentFluxes)
    print(io, summary(fluxes), '\n',
          "├── gravitational_acceleration: ", prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",        prettysummary(fluxes.von_karman_constant), '\n',
          "├── turbulent_prandtl_number: ",   prettysummary(fluxes.turbulent_prandtl_number), '\n',
          "├── gustiness_parameter: ",        prettysummary(fluxes.gustiness_parameter), '\n',
          "├── stability_functions: ",        summary(fluxes.stability_functions), '\n',
          "├── water_mole_fraction: ",        summary(fluxes.water_mole_fraction), '\n',
          "├── water_vapor_saturation: ",     summary(fluxes.water_vapor_saturation), '\n',
          "├── roughness_lengths: ",          summary(fluxes.roughness_lengths), '\n',
          "├── similarity_profile_type: ",    summary(fluxes.similarity_profile_type), '\n',
          "├── surface_temperature: ",        summary(fluxes.surface_temperature_type), '\n',
          "└── thermodynamics_parameters: ",  summary(fluxes.thermodynamics_parameters))
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

""" The exchange fluxes depend on the atmosphere velocity but not the ocean velocity """
struct WindVelocity end

""" The exchange fluxes depend on the relative velocity between the atmosphere and the ocean """
struct RelativeVelocity end

"""
    SimilarityTheoryTurbulentFluxes(FT::DataType = Float64;
                                    gravitational_acceleration = default_gravitational_acceleration,
                                    von_karman_constant = convert(FT, 0.4),
                                    turbulent_prandtl_number = convert(FT, 1),
                                    gustiness_parameter = convert(FT, 6.5),
                                    stability_functions = default_stability_functions(FT),
                                    thermodynamics_parameters = PATP(FT),
                                    water_vapor_saturation = ClasiusClapyeronSaturation(),
                                    water_mole_fraction = convert(FT, 0.98),
                                    roughness_lengths = default_roughness_lengths(FT),
                                    similarity_profile_type = LogarithmicSimilarityProfile(),
                                    surface_temperature_type = BulkTemperature(),
                                    bulk_velocity = RelativeVelocity(),
                                    tolerance = 1e-8,
                                    maxiter = 100,
                                    fields = nothing)

`SimilarityTheoryTurbulentFluxes` contains parameters and settings to calculate
sea-air turbulent fluxes using Monin-Obukhov similarity theory.

Keyword Arguments
==================

- `gravitational_acceleration`: The gravitational acceleration. Default: `default_gravitational_acceleration`.
- `von_karman_constant`: The von Karman constant. Default: 0.4.
- `turbulent_prandtl_number`: The turbulent Prandtl number. Default: 1.
- `gustiness_parameter`: The gustiness parameter that accounts for low wind speed areas. Default: 6.5.
- `stability_functions`: The stability functions. Default: `default_stability_functions(FT)` that follow the 
                         formulation of Edson et al. (2013).
- `thermodynamics_parameters`: The thermodynamics parameters used to calculate atmospheric stability and
                               saturation pressure. Default: `PATP(FT)`, alias for `PrescribedAtmosphereThermodynamicsParameters`.
- `water_vapor_saturation`: The water vapor saturation law. Default: `ClasiusClapyeronSaturation()` that follows the 
                            Clasius-Clapyeron pressure formulation.
- `water_mole_fraction`: The water mole fraction used to calculate the `seawater_saturation_specific_humidity`. 
                         Default: 0.98, the rest is assumed to be other substances such as chlorine, sodium sulfide, and magnesium.
- `roughness_lengths`: The roughness lengths used to calculate the characteristic scales for momentum, temperature and 
                       water vapor. Default: `default_roughness_lengths(FT)`, formulation taken from Edson et al (2013).
- `similarity_profile_type`: The type of similarity profile used to relate the atmospheric state to the 
                             surface fluxes / characteristic scales.
- `bulk_velocity`: The velocity used to calculate the characteristic scales. Default: `RelativeVelocity()` (difference between
                   atmospheric and oceanic speed).
- `tolerance`: The tolerance for convergence. Default: 1e-8.
- `maxiter`: The maximum number of iterations. Default: 100.
- `fields`: The fields to calculate. Default: nothing.
"""
function SimilarityTheoryTurbulentFluxes(FT::DataType = Float64;
                                         gravitational_acceleration = default_gravitational_acceleration,
                                         von_karman_constant = convert(FT, 0.4),
                                         turbulent_prandtl_number = convert(FT, 1),
                                         gustiness_parameter = convert(FT, 6.5),
                                         stability_functions = edson_stability_functions(FT),
                                         thermodynamics_parameters = PATP(FT),
                                         water_vapor_saturation = ClasiusClapyeronSaturation(),
                                         water_mole_fraction = convert(FT, 0.98),
                                         roughness_lengths = default_roughness_lengths(FT),
                                         similarity_profile_type = LogarithmicSimilarityProfile(),
                                         surface_temperature_type = BulkTemperature(),
                                         bulk_velocity = RelativeVelocity(),
                                         tolerance = 1e-8,
                                         maxiter = 100,
                                         fields = nothing)

    return SimilarityTheoryTurbulentFluxes(convert(FT, gravitational_acceleration),
                                           convert(FT, von_karman_constant),
                                           convert(FT, turbulent_prandtl_number),
                                           convert(FT, gustiness_parameter),
                                           stability_functions,
                                           thermodynamics_parameters,
                                           water_vapor_saturation,
                                           water_mole_fraction,
                                           roughness_lengths,
                                           similarity_profile_type,
                                           surface_temperature_type,
                                           bulk_velocity,
                                           convert(FT, tolerance), 
                                           maxiter,
                                           fields)
end

function SimilarityTheoryTurbulentFluxes(grid::AbstractGrid; surface_temperature_type = BulkTemperature(), kw...)
    water_vapor   = Field{Center, Center, Nothing}(grid)
    latent_heat   = Field{Center, Center, Nothing}(grid)
    sensible_heat = Field{Center, Center, Nothing}(grid)
    x_momentum    = Field{Center, Center, Nothing}(grid)
    y_momentum    = Field{Center, Center, Nothing}(grid)
    T_surface     = Field{Center, Center, Nothing}(grid)

    fields = (; latent_heat, sensible_heat, water_vapor, x_momentum, y_momentum, T_surface)

    return SimilarityTheoryTurbulentFluxes(eltype(grid); surface_temperature_type, kw..., fields)
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

#####
##### Fixed-point iteration for roughness length
#####

@inline function compute_similarity_theory_fluxes(similarity_theory,
                                                  surface_state,
                                                  atmos_state,
                                                  prescribed_heat_fluxes, # Possibly use in state_differences
                                                  radiative_properties,
                                                  ocean_salinity,
                                                  ocean_density,
                                                  ocean_heat_capacity,
                                                  atmos_boundary_layer_height,
                                                  thermodynamics_parameters,
                                                  gravitational_acceleration,
                                                  von_karman_constant,
                                                  maxiter)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    FT = eltype(ℂₐ)

    # Initial guess for the characteristic scales u★, θ★, q★.
    # Does not really matter if we are sophisticated or not, it converges 
    # in about 10 iterations no matter what...
    Δu, Δv = velocity_differences(atmos_state, surface_state, similarity_theory.bulk_velocity)

    # The inital velocity scale assumes that the gustiness velocity `Uᴳ` is equal to 0.5 ms⁻¹. 
    # The initial surface temperature is the same as the ocean temperature.
    # These will be refined later on.
    θs   = AtmosphericThermodynamics.air_temperature(ℂₐ, surface_state.ts)
    Uᴳᵢ² = convert(FT, 0.5^2)
    ΔU   = sqrt(Δu^2 + Δv^2 + Uᴳᵢ²)
    
    # break the cycle if Δu == Δv == gustiness_parameter == 0 since that would mean 
    # that u★ == 0 so there is no turbulent transfer and the solver will not converge, leading to NaNs.
    zero_shear_velocity = (Δu == 0) & (Δv == 0) & (similarity_theory.gustiness_parameter == 0)

    # Initialize the solver
    iteration = ifelse(zero_shear_velocity, maxiter+1, 0)
    u★ = ifelse(zero_shear_velocity, zero(FT), convert(FT, 1e-4))
    Σ★ = SimilarityScales(u★, u★, u★) 
    Σ₀ = Σ★

    # Iterate until convergence
    while iterating(Σ★ - Σ₀, iteration, maxiter, similarity_theory)
        Σ₀ = Σ★
        # Refine both the characteristic scale, the effective
        # velocity difference ΔU, including gustiness, and the surface
        # state temperature.
        Σ★, θs, ΔU = refine_similarity_variables(Σ★, θs, ΔU,
                                                 similarity_theory,
                                                 atmos_state,
                                                 surface_state,
                                                 ocean_salinity,
                                                 ocean_density,
                                                 ocean_heat_capacity,
                                                 atmos_boundary_layer_height,
                                                 thermodynamics_parameters,
                                                 prescribed_heat_fluxes,
                                                 radiative_properties,
                                                 gravitational_acceleration,
                                                 von_karman_constant)
        iteration += 1
    end

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

    θ★ = θ★ / similarity_theory.turbulent_prandtl_number
    q★ = q★ / similarity_theory.turbulent_prandtl_number

    # `u★² ≡ sqrt(τx² + τy²)`
    # We remove the gustiness by dividing by `ΔU`
    τx = - u★^2 * Δu / ΔU
    τy = - u★^2 * Δv / ΔU

    𝒬ₐ = atmos_state.ts
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    fluxes = (;
        sensible_heat = - ρₐ * cₚ * u★ * θ★,
        latent_heat   = - ρₐ * u★ * q★ * ℰv,
        water_vapor   = - ρₐ * u★ * q★,
        x_momentum    = + ρₐ * τx,
        y_momentum    = + ρₐ * τy,
    )
    
    return fluxes, θs
end

# Iterating condition for the characteristic scales solvers
@inline function iterating(Σ★, iteration, maxiter, solver)
    hasnt_started = iteration == 0
    converged = norm(Σ★) < solver.tolerance
    reached_maxiter = iteration ≥ maxiter
    return !(converged | reached_maxiter) | hasnt_started
end

"""
    buoyancy_scale(θ★, q★, 𝒬, ℂ, g)

Return the characteristic buoyancy scale `b★` associated with
the characteristic temperature `θ★`, specific humidity scale `q★`,
near-surface atmospheric thermodynamic state `𝒬`, thermodynamic
parameters `ℂ`, and gravitational acceleration `g`.

The buoyancy scale is defined in terms of the surface buoyancy flux,

```math
u★ b★ ≡ w′b′,
```

where `u★` is the friction velocity.
Using the definition of buoyancy for non-condensing air, we find that

```math
b★ = g / 𝒯ₐ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₐ * q★),
```
where ``𝒯ₐ`` is the virtual temperature of the atmosphere near the surface,
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

@inline velocity_differences(𝒰₁, 𝒰₀, ::RelativeVelocity) = @inbounds 𝒰₁.u[1] - 𝒰₀.u[1], 𝒰₁.u[2] - 𝒰₀.u[2]
@inline velocity_differences(𝒰₁, 𝒰₀, ::WindVelocity)     = @inbounds 𝒰₁.u[1], 𝒰₁.u[2] 

@inline function state_differences(ℂ, 𝒰₁, 𝒰₀, θ₀, S₀, Σ★, g, ρₒ, cpₒ, 
                                   water_mole_fraction,
                                   water_vapor_saturation,
                                   surface_temperature_type, 
                                   prescribed_heat_fluxes,
                                   radiative_properties,
                                   bulk_velocity)
    z₁ = 𝒰₁.z
    z₀ = 𝒰₀.z
    Δh = z₁ - z₀
    Δu, Δv = velocity_differences(𝒰₁, 𝒰₀, bulk_velocity)
    
    # Thermodynamic state
    𝒬₁ = 𝒰₁.ts
    𝒬₀ = 𝒰₀.ts

    ρₐ = AtmosphericThermodynamics.air_density(ℂ, 𝒬₁)
    cₚ = AtmosphericThermodynamics.cp_m(ℂ, 𝒬₁) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂ, 𝒬₁)

    θ₀ = compute_surface_temperature(surface_temperature_type, θ₀, ℂ, 𝒬₀, ρₐ, cₚ, ℰv, Σ★, ρₒ, cpₒ, g,
                                     prescribed_heat_fluxes, 
                                     radiative_properties)

    θ₁ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₁)

    # Temperature difference including the ``lapse rate'' `α = g / cₚ`
    Δθ = θ₁ - θ₀ + g / cₚ * Δh

    q₁ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₁)

    # Recomputing the saturation specific humidity at the surface based on the new temperature
    q₀ = seawater_saturation_specific_humidity(ℂ, θ₀, S₀, 𝒬₁,
                                               water_mole_fraction,
                                               water_vapor_saturation,
                                               AtmosphericThermodynamics.Liquid())
    
    𝒬ₛ = AtmosphericThermodynamics.PhaseEquil_pTq(ℂ, 𝒬₀.p, θ₀, q₀)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬ₛ)

    Δq = q₁ - q₀

    return Δh, Δu, Δv, Δθ, Δq, θ₀
end

@inline function refine_similarity_variables(estimated_characteristic_scales, 
                                             surface_temperature,
                                             velocity_scale,
                                             similarity_theory,
                                             atmos_state,
                                             surface_state,
                                             surface_salinity,
                                             ocean_density,
                                             ocean_heat_capacity,
                                             atmos_boundary_layer_height,
                                             thermodynamics_parameters,
                                             prescribed_heat_fluxes,
                                             radiative_properties,
                                             gravitational_acceleration,
                                             von_karman_constant)

    Δh, Δu, Δv, Δθ, Δq, θ₀ = state_differences(thermodynamics_parameters,
                                               atmos_state,
                                               surface_state,
                                               surface_temperature,
                                               surface_salinity,
                                               estimated_characteristic_scales,
                                               gravitational_acceleration,
                                               ocean_density,
                                               ocean_heat_capacity,
                                               similarity_theory.water_mole_fraction,
                                               similarity_theory.water_vapor_saturation,
                                               similarity_theory.surface_temperature_type,
                                               prescribed_heat_fluxes,
                                               radiative_properties,
                                               similarity_theory.bulk_velocity)
                                               
    # "initial" scales because we will recompute them
    u★ = estimated_characteristic_scales.momentum
    θ★ = estimated_characteristic_scales.temperature
    q★ = estimated_characteristic_scales.water_vapor
    ΔU = velocity_scale

    # Similarity functions from Edson et al. (2013)
    ψu = similarity_theory.stability_functions.momentum
    ψθ = similarity_theory.stability_functions.temperature
    ψq = similarity_theory.stability_functions.water_vapor

    # Extract roughness lengths
    ℓu = similarity_theory.roughness_lengths.momentum
    ℓθ = similarity_theory.roughness_lengths.temperature
    ℓq = similarity_theory.roughness_lengths.water_vapor
    β  = similarity_theory.gustiness_parameter

    ℂ  = thermodynamics_parameters
    g  = gravitational_acceleration
    𝒬ₒ = surface_state.ts # thermodynamic state

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    b★ = buoyancy_scale(θ★, q★, 𝒬ₒ, ℂ, g)

    # Monin-Obhukov characteristic length scale and non-dimensional height
    ϰ  = von_karman_constant
    L★ = ifelse(b★ == 0, zero(b★), - u★^2 / (ϰ * b★))

    # Compute roughness length scales
    ℓu₀ = roughness_length(ℓu, u★, 𝒬ₒ, ℂ)
    ℓq₀ = roughness_length(ℓq, ℓu₀, u★, 𝒬ₒ, ℂ)
    ℓθ₀ = roughness_length(ℓθ, ℓu₀, u★, 𝒬ₒ, ℂ)

    # Transfer coefficients at height `h`
    profile_type = similarity_theory.similarity_profile_type
    χu = ϰ / similarity_profile(profile_type, ψu, Δh, ℓu₀, L★)
    χθ = ϰ / similarity_profile(profile_type, ψθ, Δh, ℓθ₀, L★)
    χq = ϰ / similarity_profile(profile_type, ψq, Δh, ℓq₀, L★)

    # u★ including gustiness
    u★ = χu * ΔU
    θ★ = χθ * Δθ
    q★ = χq * Δq

    # Buoyancy flux characteristic scale for gustiness (Edson 2013)
    hᵢ = atmos_boundary_layer_height
    Jᵇ = - u★ * b★
    Uᴳ = β * cbrt(Jᵇ * hᵢ)

    # New velocity difference accounting for gustiness
    ΔU = sqrt(Δu^2 + Δv^2 + Uᴳ^2)

    return SimilarityScales(u★, θ★, q★), θ₀, ΔU
end
