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
    universal_func_type,
    grav


#####
##### Bulk turbulent fluxes based on similarity theory
#####

struct SimilarityTheoryTurbulentFluxes{FT, UF, TP, S, W, R, B, V, F}
    gravitational_acceleration :: FT
    von_karman_constant :: FT
    turbulent_prandtl_number :: FT
    gustiness_parameter :: FT
    stability_functions :: UF
    thermodynamics_parameters :: TP
    water_vapor_saturation :: S
    water_mole_fraction :: W
    roughness_lengths :: R
    bulk_coefficients :: B
    bulk_velocity :: V
    tolerance :: FT
    maxiter :: Int
    fields :: F
end

const STTF = SimilarityTheoryTurbulentFluxes
@inline thermodynamics_params(fluxes::STTF) = fluxes.thermodynamics_parameters
@inline uf_params(fluxes::STTF)             = fluxes.stability_functions
@inline von_karman_const(fluxes::STTF)      = fluxes.von_karman_constant
@inline grav(fluxes::STTF)                  = fluxes.gravitational_acceleration
@inline molmass_ratio(fluxes::STTF)         = molmass_ratio(fluxes.thermodynamics_parameters)

@inline universal_func_type(::STTF{<:Any, <:Any, <:BusingerParams}) = BusingerType()

Adapt.adapt_structure(to, fluxes::STTF) = SimilarityTheoryTurbulentFluxes(adapt(to, fluxes.gravitational_acceleration),
                                                                          adapt(to, fluxes.von_karman_constant),
                                                                          adapt(to, fluxes.turbulent_prandtl_number),
                                                                          adapt(to, fluxes.gustiness_parameter),
                                                                          adapt(to, fluxes.stability_functions),
                                                                          adapt(to, fluxes.thermodynamics_parameters),
                                                                          adapt(to, fluxes.water_vapor_saturation),
                                                                          adapt(to, fluxes.water_mole_fraction),
                                                                          adapt(to, fluxes.roughness_lengths),
                                                                          adapt(to, fluxes.bulk_coefficients),
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
          "├── gravitational_acceleration: ",      prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",             prettysummary(fluxes.von_karman_constant), '\n',
          "├── turbulent_prandtl_number: ",        prettysummary(fluxes.turbulent_prandtl_number), '\n',
          "├── gustiness_parameter: ",             prettysummary(fluxes.gustiness_parameter), '\n',
          "├── stability_functions: ",             summary(fluxes.stability_functions), '\n',
          "├── water_mole_fraction: ",             summary(fluxes.water_mole_fraction), '\n',
          "├── water_vapor_saturation: ",          summary(fluxes.water_vapor_saturation), '\n',
          "├── roughness_lengths: ",               summary(fluxes.roughness_lengths), '\n',
          "├── bulk_coefficients: ",               summary(fluxes.bulk_coefficients), '\n',
          "└── thermodynamics_parameters: ",       summary(fluxes.thermodynamics_parameters))
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
                                    bulk_coefficients = bulk_coefficients,
                                    bulk_velocity = RelativeVelocity(),
                                    tolerance = 1e-8,
                                    maxiter = 100,
                                    fields = nothing)

`SimilarityTheoryTurbulentFluxes` contains parameters and settings to calculate
sea-air turbulent fluxes using Monin-Obukhov similarity theory.

Keyword Arguments
==================

- `gravitational_acceleration`: The gravitational acceleration (default: default_gravitational_acceleration).
- `von_karman_constant`: The von Karman constant (default: 0.4).
- `turbulent_prandtl_number`: The turbulent Prandtl number (default: 1).
- `gustiness_parameter`: The gustiness parameter that accounts for low wind speed areas (default: 6.5).
- `stability_functions`: The stability functions. Default: default_stability_functions(FT) that follow the 
                         formulation of Edson et al (2013).
- `thermodynamics_parameters`: The thermodynamics parameters used to calculate atmospheric stability and
                               saturation pressure. Default: `PATP`, alias for `PrescribedAtmosphereThermodynamicsParameters`.
- `water_vapor_saturation`: The water vapor saturation law. Default: ClasiusClapyeronSaturation() that follows the 
                            Clasius Clapyeron pressure formulation.
- `water_mole_fraction`: The water mole fraction used to calculate the seawater_saturation_specific_humidity. 
                         Default: 0.98, the rest is assumed to be other substances such as chlorine, sodium sulfide and magnesium.
- `roughness_lengths`: The roughness lengths used to calculate the characteristic scales for momentum, temperature and 
                       water vapor. Default: default_roughness_lengths(FT), formulation taken from Edson et al (2013).
- `bulk_coefficients`: The bulk coefficients.
- `bulk_velocity`: The velocity used to calculate the characteristic scales. Default: RelativeVelocity() (difference between
                   atmospheric and oceanic speed).
- `tolerance`: The tolerance for convergence (default: 1e-8).
- `maxiter`: The maximum number of iterations (default: 100).
- `fields`: The fields to calculate (default: nothing).
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
                                         bulk_coefficients = bulk_coefficients,
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
                                           bulk_coefficients,
                                           bulk_velocity,
                                           convert(FT, tolerance), 
                                           maxiter,
                                           fields)
end

function SimilarityTheoryTurbulentFluxes(grid::AbstractGrid; kw...)
    water_vapor   = Field{Center, Center, Nothing}(grid)
    latent_heat   = Field{Center, Center, Nothing}(grid)
    sensible_heat = Field{Center, Center, Nothing}(grid)
    x_momentum    = Field{Center, Center, Nothing}(grid)
    y_momentum    = Field{Center, Center, Nothing}(grid)

    fields = (; latent_heat, sensible_heat, water_vapor, x_momentum, y_momentum)

    return SimilarityTheoryTurbulentFluxes(eltype(grid); kw..., fields)
end

# Simplified coefficient a la COARE 
@inline simplified_bulk_coefficients(ψ, h, ℓ, L) = log(h / ℓ) - ψ(h / L) # + ψ(ℓ / L)

# The complete bulk coefficient
@inline bulk_coefficients(ψ, h, ℓ, L) = log(h / ℓ) - ψ(h / L) + ψ(ℓ / L)

#####
##### Fixed-point iteration for roughness length
#####

@inline function compute_similarity_theory_fluxes(similarity_theory,
                                                  surface_state,
                                                  atmos_state,
                                                  atmos_boundary_layer_height,
                                                  thermodynamics_parameters,
                                                  gravitational_acceleration,
                                                  von_karman_constant,
                                                  maxiter)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    Δh, Δu, Δv, Δθ, Δq = state_differences(ℂₐ, 
                                           atmos_state, 
                                           surface_state, 
                                           gravitational_acceleration,
                                           similarity_theory.bulk_velocity)

    differences = (; u=Δu, v=Δv, θ=Δθ, q=Δq, h=Δh)
    
    u★ = convert(eltype(Δh), 1e-4)

    # Initial guess for the characteristic scales u★, θ★, q★.
    # Does not really matter if we are sophisticated or not, it converges 
    # in about 10 iterations no matter what...
    Σ₀ = SimilarityScales(1, 1, 1)
    Σ★ = SimilarityScales(u★, u★, u★) 

    # The inital velocity scale assumes that
    # the gustiness velocity `uᴳ` is equal to 0.5 ms⁻¹. 
    # That will be refined later on.
    uτ = sqrt(Δu^2 + Δv^2 + convert(eltype(Δh), 0.25))

    # Initialize the solver
    iteration = 0

    while iterating(Σ★ - Σ₀, iteration, maxiter, similarity_theory)
        Σ₀ = Σ★
        Σ★, uτ, = refine_characteristic_scales(Σ★, uτ, 
                                               similarity_theory,
                                               surface_state,
                                               differences,
                                               atmos_boundary_layer_height,
                                               thermodynamics_parameters,
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
    # We remove the gustiness by dividing by `uτ`
    τx = - u★^2 * Δu / uτ
    τy = - u★^2 * Δv / uτ

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

    return fluxes
end

# Iterating condition for the characteristic scales solvers
@inline function iterating(Σ★, iteration, maxiter, solver)
    converged = norm(Σ★) <= solver.tolerance
    reached_maxiter = iteration >= maxiter 
    return !(converged | reached_maxiter)
end

# The M-O characteristic length is calculated as
#  L★ = - u★² / (κ ⋅ b★)
# where b★ is the characteristic buoyancy scale calculated from:
@inline function buoyancy_scale(θ★, q★, 𝒬, ℂ, g)
    𝒯ₐ = AtmosphericThermodynamics.virtual_temperature(ℂ, 𝒬)
    qₐ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬)
    ε  = AtmosphericThermodynamics.Parameters.molmass_ratio(ℂ)
    δ  = ε - 1 # typically equal to 0.608

    # Fairell et al. 1996, 
    b★ = g / 𝒯ₐ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₐ * q★)

    return b★
end

@inline velocity_differences(𝒰₁, 𝒰₀, ::RelativeVelocity) = @inbounds 𝒰₁.u[1] - 𝒰₀.u[1], 𝒰₁.u[2] - 𝒰₀.u[2]
@inline velocity_differences(𝒰₁, 𝒰₀, ::WindVelocity)     = @inbounds 𝒰₁.u[1], 𝒰₁.u[2] 

@inline function state_differences(ℂ, 𝒰₁, 𝒰₀, g, bulk_velocity)
    z₁ = 𝒰₁.z
    z₀ = 𝒰₀.z
    Δh = z₁ - z₀
    Δu, Δv = velocity_differences(𝒰₁, 𝒰₀, bulk_velocity)

    # Thermodynamic state
    𝒬₁ = 𝒰₁.ts
    𝒬₀ = 𝒰₀.ts

    θ₁ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₁)
    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)
    cₚ = AtmosphericThermodynamics.cp_m(ℂ, 𝒬₁) # moist heat capacity

    # Temperature difference including the ``lapse rate'' `α = g / cₚ`
    Δθ = θ₁ - θ₀ + g / cₚ * Δh

    q₁ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₁)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₀)
    Δq = q₁ - q₀

    return Δh, Δu, Δv, Δθ, Δq
end

@inline function refine_characteristic_scales(estimated_characteristic_scales, 
                                              velocity_scale,
                                              similarity_theory,
                                              surface_state,
                                              differences,
                                              atmos_boundary_layer_height,
                                              thermodynamics_parameters,
                                              gravitational_acceleration,
                                              von_karman_constant)

    # "initial" scales because we will recompute them
    u★ = estimated_characteristic_scales.momentum
    θ★ = estimated_characteristic_scales.temperature
    q★ = estimated_characteristic_scales.water_vapor
    uτ = velocity_scale

    # Similarity functions from Edson et al. (2013)
    ψu = similarity_theory.stability_functions.momentum
    ψθ = similarity_theory.stability_functions.temperature
    ψq = similarity_theory.stability_functions.water_vapor

    # Extract roughness lengths
    ℓu = similarity_theory.roughness_lengths.momentum
    ℓθ = similarity_theory.roughness_lengths.temperature
    ℓq = similarity_theory.roughness_lengths.water_vapor
    β  = similarity_theory.gustiness_parameter

    h  = differences.h
    ϰ  = von_karman_constant
    ℂ  = thermodynamics_parameters
    g  = gravitational_acceleration
    𝒬ₒ = surface_state.ts # thermodynamic state
    zᵢ = atmos_boundary_layer_height

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    b★ = buoyancy_scale(θ★, q★, 𝒬ₒ, ℂ, g)

    # Monin-Obhukov characteristic length scale and non-dimensional height
    L★ = ifelse(b★ == 0, zero(b★), - u★^2 / (ϰ * b★))
    
    # Compute roughness length scales
    ℓu₀ = roughness_length(ℓu, u★, 𝒬ₒ, ℂ)
    ℓq₀ = roughness_length(ℓq, ℓu₀, u★, 𝒬ₒ, ℂ)
    ℓθ₀ = roughness_length(ℓθ, ℓu₀, u★, 𝒬ₒ, ℂ)

    # Transfer coefficients at height `h`
    χu = ϰ / similarity_theory.bulk_coefficients(ψu, h, ℓu₀, L★) 
    χθ = ϰ / similarity_theory.bulk_coefficients(ψθ, h, ℓθ₀, L★) 
    χq = ϰ / similarity_theory.bulk_coefficients(ψq, h, ℓq₀, L★) 

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    # u★ including gustiness
    u★ = χu * uτ
    θ★ = χθ * Δθ
    q★ = χq * Δq

    # Buoyancy flux characteristic scale for gustiness (Edson 2013)
    Jᵇ = - u★ * b★
    uᴳ = β * cbrt(Jᵇ * zᵢ)

    # New velocity difference accounting for gustiness
    uτ = sqrt(Δu^2 + Δv^2 + uᴳ^2)

    return SimilarityScales(u★, θ★, q★), uτ
end