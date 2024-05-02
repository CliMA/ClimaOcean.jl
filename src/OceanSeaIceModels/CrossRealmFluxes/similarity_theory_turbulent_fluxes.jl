using Oceananigans.Utils: prettysummary
using Oceananigans.Grids: AbstractGrid

using Adapt
using Thermodynamics: Liquid
using SurfaceFluxes.Parameters: SurfaceFluxesParameters, AbstractSurfaceFluxesParameters
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

struct SimilarityTheoryTurbulentFluxes{FT, UF, TP, S, W, R, F} <: AbstractSurfaceFluxesParameters
    gravitational_acceleration :: FT
    von_karman_constant :: FT
    turbulent_prandtl_number :: FT
    gustiness_parameter :: FT
    stability_functions :: UF
    thermodynamics_parameters :: TP
    water_vapor_saturation :: S
    water_mole_fraction :: W
    roughness_lengths :: R
    tolerance :: FT
    maxiter :: Int
    iteration :: Int
    fields :: F
end

const STTF = SimilarityTheoryTurbulentFluxes
@inline thermodynamics_params(fluxes::STTF) = fluxes.thermodynamics_parameters
@inline uf_params(fluxes::STTF)             = fluxes.stability_functions
@inline von_karman_const(fluxes::STTF)      = fluxes.von_karman_constant
@inline grav(fluxes::STTF)                  = fluxes.gravitational_acceleration
@inline molmass_ratio(fluxes::STTF)         = molmass_ratio(fluxes.thermodynamics_parameters)

@inline universal_func_type(fluxes::STTF{<:Any, <:Any, <:BusingerParams}) = BusingerType()

Adapt.adapt_structure(to, fluxes::STTF) = SimilarityTheoryTurbulentFluxes(adapt(to, fluxes.gravitational_acceleration),
                                                                          adapt(to, fluxes.von_karman_constant),
                                                                          adapt(to, fluxes.turbulent_prandtl_number),
                                                                          adapt(to, fluxes.gustiness_parameter),
                                                                          adapt(to, fluxes.stability_functions),
                                                                          adapt(to, fluxes.thermodynamics_parameters),
                                                                          adapt(to, fluxes.water_vapor_saturation),
                                                                          adapt(to, fluxes.water_mole_fraction),
                                                                          adapt(to, fluxes.roughness_lengths),
                                                                          fluxes.iteration,
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
          "└── thermodynamics_parameters: ",       summary(fluxes.thermodynamics_parameters))
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

function SimilarityTheoryTurbulentFluxes(FT::DataType = Float64;
                                         gravitational_acceleration = default_gravitational_acceleration,
                                         von_karman_constant = convert(FT, 0.4),
                                         turbulent_prandtl_number = convert(FT, 1),
                                         gustiness_parameter = convert(FT, 6.5),
                                         stability_functions = default_stability_functions(FT),
                                         thermodynamics_parameters = PATP(FT),
                                         water_vapor_saturation = ClasiusClapyeronSaturation(),
                                         water_mole_fraction = convert(FT, 0.98),
                                         roughness_lengths = default_roughness_lengths(FT),
                                         tolerance = 1e-8,
                                         maxiter = Inf,
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
                                           tolerance, 
                                           maxiter,
                                           0,
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

#####
##### Fixed-point iteration for roughness length
#####

@inline function compute_similarity_theory_fluxes(similarity_theory,
                                                  surface_state,
                                                  atmos_state,
                                                  atmos_boundary_layer_height,
                                                  thermodynamics_parameters,
                                                  gravitational_acceleration,
                                                  von_karman_constant)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    Δh, Δu, Δv, Δθ, Δq = state_differences(ℂₐ, atmos_state, surface_state, gravitational_acceleration)
    differences = (; u=Δu, v=Δv, θ=Δθ, q=Δq, h=Δh)
    
    Σ₀ = SimilarityScales(0, 0, 0)

    # Initial guess for the characteristic scales u★, θ★, q★.
    Σ★ = initial_guess(differences, 
                       similarity_theory,
                       atmos_boundary_layer_height,
                       gravitational_acceleration,
                       von_karman_constant, 
                       ℂₐ, 
                       surface_state.ts)

    # The inital velocity scale assumes that
    # the gustiness velocity `uᴳ` is equal to 0.5 ms⁻¹. 
    # That will be refined later on.
    uτ = sqrt(Δu^2 + Δv^2 + convert(FT, 0.25))

    while iterating(Σ★ - Σ₀, similarity_theory)
        Σ₀ = Σ★
        similarity_theory.iteration += 1
        Σ★, uτ, = refine_characteristic_scales(Σ★, uτ, 
                                               similarity_theory,
                                               surface_state,
                                               differences,
                                               atmos_boundary_layer_height,
                                               thermodynamics_parameters,
                                               gravitational_acceleration,
                                               von_karman_constant)
    end

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

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
@inline function iterating(Σ★, solver)
    solver.iteration >= solver.maxiter && return false
    norm(Σ★) <= solver.tolerance && return false
    return true
end

# The complete bulk coefficient should include also `ψ(ℓ / L)`, but the 
# JRA55 atmosphere is adjusted to formulae without this last term so we exclude it
@inline bulk_coefficient(ψ, h, ℓ, L) = log(h / ℓ) - ψ(h / L) # + ψ(ℓ / L)

@inline function initial_guess(differences, 
                               similarity_theory,
                               atmos_boundary_layer_height,
                               gravitational_acceleration,
                               von_karman_constant, 
                               ℂₐ, 𝒬ₒ)

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q
    h  = differences.h

    g  = gravitational_acceleration
    ϰ  = von_karman_constant

    # Extract roughness lengths
    ℓu  = similarity_theory.roughness_lengths.momentum
    β   = similarity_theory.gustiness_parameter
    zᵢ  = atmos_boundary_layer_height

    hᵢ  = convert(eltype(h), 10)    # Reference Initial height == 10 meters
    ℓuᵢ = convert(eltype(h), 1e-4)  # Initial roughness length == 1e-4 meters

    # assuming the initial gustiness is `0.5` ms⁻¹
    uτ = sqrt(Δu^2 + Δv^2 + convert(FT, 0.25))

    # u10 at the reference ten meter height, assuming the initial roughness length is `1e-4` m
    u10 = uτ / log(h / ℓuᵢ) * 11.5129 # log(10 / 1e-4) == 11.5129
    u★  = convert(FT, 0.035) * u10

    ℓu₀ = roughness_length(ℓu, u★, 𝒬ₐ, ℂₐ)

    # Initial neutral coefficients at 10 meter height
    χuₙ  = (ϰ / log(hᵢ / ℓu₀))^2
    χcₙ  = convert(FT, 0.00115) / sqrt(χuₙ)

    # Initial scalar roughness length
    ℓθ₀ = hᵢ / exp(ϰ / χcₙ)

    # Neutral transfer coefficients at height `h`
    χu = (ϰ / log(h / ℓu₀))^2
    χq =  ϰ / log(h / ℓθ₀)
    χc =  ϰ * χq / χu
    
    # Similarity functions from Businger et al. (1971)
    ψu = InitialMomentumStabilityFunction()
    ψθ = similarity_theory.stability_functions.temperature
    ψq = similarity_theory.stability_functions.water_vapor

    # Bulk Flux Richardson number
    # TODO: find out what 0.004 refers to
    b★  = buoyancy_scale(Δθ, Δq, 𝒬ₒ, ℂₐ, g)
    Ri  = - ifelse(b★ == 0, zero(b★), h / b★ / uτ^2)
    Riᶜ = - h / zᵢ / convert(FT, 0.004) / β^3 # - h / zi / 0.004 / β^3
    
    # Calculating the first stability coefficient and the MO length
    # TODO: explain this formulation of the stability function. Is it empirical?
    # Found in COARE3.6
    ζ10 = ifelse(Ri < 0, χc * Ri / (1 + Ri / Riᶜ), χc * Ri * (1 + 27 / 9 * Ri / χc))
    L10 = h / ζ10

    u★ = uτ * ϰ / bulk_coefficient(ψu, h, ℓu₀, L10) 
    θ★ = Δθ * ϰ / bulk_coefficient(ψθ, h, ℓθ₀, L10) 
    q★ = Δq * ϰ / bulk_coefficient(ψq, h, ℓθ₀, L10) 
    
    return SimilarityScales(u★, θ★, q★)
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

@inline function state_differences(ℂ, 𝒰₁, 𝒰₀, g)
    z₁ = 𝒰₁.z
    z₀ = 𝒰₀.z
    Δh = z₁ - z₀

    U₁ = 𝒰₁.u
    U₀ = 𝒰₀.u

    @inbounds begin
        Δu = U₁[1] - U₀[1]
        Δv = U₁[2] - U₀[2]
    end

    # Thermodynamic state
    𝒬₁ = 𝒰₁.ts
    𝒬₀ = 𝒰₀.ts

    θ₁ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₁)
    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₀)
    cₚ = AtmosphericThermodynamics.cp_m(ℂ, 𝒬₁) # moist heat capacity

    # The temperature difference includes the ``lapse rate'' α = g / h
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

    # Similarity functions from Edison et al. (2013)
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
    χu = ϰ / bulk_coefficient(ψu, h, ℓu₀, L★) 
    χθ = ϰ / bulk_coefficient(ψθ, h, ℓθ₀, L★) 
    χq = ϰ / bulk_coefficient(ψq, h, ℓq₀, L★) 

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    # u★ including gustiness
    u★ = χu * uτ
    θ★ = χθ * Δθ
    q★ = χq * Δq

    # Buoyancy flux characteristic scale for gustiness (Edison 2013)
    ε★ = - u★ * b★
    uᴳ = β * cbrt(ε★ * zᵢ)

    # New velocity difference accounting for gustiness
    uτ = sqrt(Δu^2 + Δv^2 + uᴳ^2)

    return SimilarityScales(u★, θ★, q★), uτ
end
