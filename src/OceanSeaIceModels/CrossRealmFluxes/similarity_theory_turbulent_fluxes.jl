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

struct SimilarityTheoryFluxes{FT, UF, R, B, T, V}
    von_karman_constant :: FT        # parameter
    turbulent_prandtl_number :: FT   # parameter
    gustiness_parameter :: FT        # bulk velocity parameter
    stability_functions :: UF        # functions for turbulent fluxes
    roughness_lengths :: R           # parameterization for turbulent fluxes
    similarity_profile_type :: B     # similarity profile relating atmosphere to surface state
    surface_temperature_type :: T    # surface temperature either diagnostic or prescribed
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
                           adapt(to, fluxes.surface_temperature_type),
                           adapt(to, fluxes.bulk_velocity),
                           fluxes.tolerance,
                           fluxes.maxiter)

Base.summary(::SimilarityTheoryFluxes{FT}) where FT = "SimilarityTheoryFluxes{$FT}"

function Base.show(io::IO, fluxes::SimilarityTheoryFluxes)
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

"""
    SimilarityTheoryFluxes(FT::DataType = Float64;
                           gravitational_acceleration = default_gravitational_acceleration,
                           von_karman_constant = convert(FT, 0.4),
                           turbulent_prandtl_number = convert(FT, 1),
                           gustiness_parameter = convert(FT, 6.5),
                           stability_functions = default_stability_functions(FT),
                           roughness_lengths = default_roughness_lengths(FT),
                           similarity_profile_type = LogarithmicSimilarityProfile(),
                           surface_temperature_type = BulkTemperature(),
                           bulk_velocity = RelativeVelocity(),
                           tolerance = 1e-8,
                           maxiter = 100)

`SimilarityTheoryFluxes` contains parameters and settings to calculate
surface-air turbulent fluxes using Monin-Obukhov similarity theory.

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
                             surface fluxes / characteristic scales.
- `bulk_velocity`: The velocity used to calculate the characteristic scales. Default: `RelativeVelocity()` (difference between
                   atmospheric and surfaceic speed).
- `tolerance`: The tolerance for convergence. Default: 1e-8.
- `maxiter`: The maximum number of iterations. Default: 100.
"""
function SimilarityTheoryFluxes(FT::DataType = Float64;
                                von_karman_constant = convert(FT, 0.4),
                                turbulent_prandtl_number = convert(FT, 1),
                                gustiness_parameter = convert(FT, 6.5),
                                stability_functions = edson_stability_functions(FT),
                                roughness_lengths = default_roughness_lengths(FT),
                                similarity_profile_type = LogarithmicSimilarityProfile(),
                                surface_temperature_type = BulkTemperature(),
                                bulk_velocity = RelativeVelocity(),
                                tolerance = 1e-8,
                                maxiter = 100)

    return SimilarityTheoryFluxes(convert(FT, von_karman_constant),
                                  convert(FT, turbulent_prandtl_number),
                                  convert(FT, gustiness_parameter),
                                  stability_functions,
                                  roughness_lengths,
                                  similarity_profile_type,
                                  surface_temperature_type,
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

#####
##### Fixed-point iteration for roughness length
#####

@inline function compute_turbulent_fluxes(similarity_theory::SimilarityTheoryFluxes,
                                          surface_state,
                                          atmos_state,
                                          prescribed_heat_fluxes, # Possibly use in state_differences
                                          radiative_properties,
                                          surface_phase,
                                          surface_salinity,
                                          surface_density,
                                          surface_heat_capacity,
                                          mole_fraction,
                                          vapor_saturation,
                                          atmos_boundary_layer_height,
                                          thermodynamics_parameters,
                                          gravitational_acceleration,
                                          inactive_cell)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    FT = eltype(ℂₐ)

    von_karman_constant = similarity_theory.von_karman_constant
    maxiter = ifelse(inactive_cell, 1, similarity_theory.maxiter)

    # Initial guess for the characteristic scales u★, θ★, q★.
    # Does not really matter if we are sophisticated or not, it converges 
    # in about 10 iterations no matter what...
    Δu, Δv = velocity_differences(atmos_state, surface_state, similarity_theory.bulk_velocity)

    # The inital velocity scale assumes that the gustiness velocity `Uᴳ` is equal to 0.5 ms⁻¹. 
    # The initial surface temperature is the same as the surface temperature.
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
                                                 surface_phase,
                                                 surface_salinity,
                                                 surface_density,
                                                 surface_heat_capacity,
                                                 mole_fraction,
                                                 vapor_saturation,
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

@inline function refine_similarity_variables(estimated_characteristic_scales, 
                                             surface_temperature,
                                             velocity_scale,
                                             similarity_theory,
                                             atmos_state,
                                             surface_state,
                                             surface_phase, # Either liquid or solid
                                             surface_salinity,
                                             surface_density,
                                             surface_heat_capacity,
                                             mole_fraction,
                                             vapor_saturation,
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
                                               surface_density,
                                               surface_heat_capacity,
                                               mole_fraction,
                                               vapor_saturation,
                                               similarity_theory.surface_temperature_type,
                                               prescribed_heat_fluxes,
                                               radiative_properties,
                                               similarity_theory.bulk_velocity,
                                               surface_phase)
                                               
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
