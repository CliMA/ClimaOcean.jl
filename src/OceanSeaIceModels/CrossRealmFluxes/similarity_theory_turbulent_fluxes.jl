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

struct SimilarityTheoryTurbulentFluxes{FT, ΔU, UF, TP, S, W, R, F} <: AbstractSurfaceFluxesParameters
    gravitational_acceleration :: FT
    von_karman_constant :: FT
    planetary_boundary_layer_height :: FT
    turbulent_prandtl_number :: FT
    bulk_velocity_scale :: ΔU
    similarity_functions :: UF
    thermodynamics_parameters :: TP
    water_vapor_saturation :: S
    water_mole_fraction :: W
    roughness_lengths :: R
    fields :: F
end

const STTF = SimilarityTheoryTurbulentFluxes
@inline thermodynamics_params(fluxes::STTF) = fluxes.thermodynamics_parameters
@inline uf_params(fluxes::STTF)             = fluxes.similarity_functions
@inline von_karman_const(fluxes::STTF)      = fluxes.von_karman_constant
@inline grav(fluxes::STTF)                  = fluxes.gravitational_acceleration
@inline molmass_ratio(fluxes::STTF)         = molmass_ratio(fluxes.thermodynamics_parameters)

@inline universal_func_type(fluxes::STTF{<:Any, <:Any, <:BusingerParams}) = BusingerType()

Adapt.adapt_structure(to, fluxes::STTF) = SimilarityTheoryTurbulentFluxes(adapt(to, fluxes.gravitational_acceleration),
                                                                          adapt(to, fluxes.von_karman_constant),
                                                                          adapt(to, fluxes.planetary_boundary_layer_height),
                                                                          adapt(to, fluxes.turbulent_prandtl_number),
                                                                          nothing, # adapt(to, fluxes.bulk_velocity_scale),
                                                                          adapt(to, fluxes.similarity_functions),
                                                                          adapt(to, fluxes.thermodynamics_parameters),
                                                                          nothing, #adapt(to, fluxes.water_vapor_saturation),
                                                                          nothing, #adapt(to, fluxes.water_mole_fraction),
                                                                          adapt(to, fluxes.roughness_lengths),
                                                                          adapt(to, fluxes.fields))

Base.summary(::SimilarityTheoryTurbulentFluxes{FT}) where FT = "SimilarityTheoryTurbulentFluxes{$FT}"

struct ClasiusClapyeronSaturation end
 
@inline function water_saturation_specific_humidity(::ClasiusClapyeronSaturation, ℂₐ, ρₛ, Tₛ)
    FT = eltype(ℂₐ)
    p★ = AtmosphericThermodynamics.saturation_vapor_pressure(ℂₐ, convert(FT, Tₛ), Liquid())
    q★ = AtmosphericThermodynamics.q_vap_saturation_from_density(ℂₐ, convert(FT, Tₛ), ρₛ, p★)
    return q★
end

struct LargeYeagerSaturation{FT}
    c₁ :: FT
    c₂ :: FT
end

function LargeYeagerSaturation(FT=Float64; c₁ = 640380, c₂ = 5107.4)
    return LargeYeagerSaturation(convert(FT, c₁), convert(FT, c₂))
end

const LYS = LargeYeagerSaturation
@inline water_saturation_specific_humidity(lys::LYS, ℂₐ, ρₛ, Tₛ) = lys.c₁ * exp(-lys.c₂ / Tₛ) / ρₛ

function Base.show(io::IO, fluxes::SimilarityTheoryTurbulentFluxes)
    print(io, summary(fluxes), '\n',
          "├── gravitational_acceleration: ",      prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",             prettysummary(fluxes.von_karman_constant), '\n',
          "├── bulk_velocity_scale: ",             summary(fluxes.bulk_velocity_scale), '\n',
          "├── planetary_boundary_layer_height: ", prettysummary(fluxes.planetary_boundary_layer_height), '\n',
          "├── turbulent_prandtl_number: ",        prettysummary(fluxes.turbulent_prandtl_number), '\n',
          "├── similarity_function: ",             summary(fluxes.similarity_function), '\n',
          "├── water_mole_fraction: ",             summary(fluxes.water_mole_fraction), '\n',
          "├── water_vapor_saturation: ",          summary(fluxes.water_vapor_saturation), '\n',
          "└── thermodynamics_parameters: ",       summary(fluxes.thermodynamics_parameters))
end

function default_roughness_lengths(FT=Float64)
    momentum    = GravityMomentumRoughnessLength(FT)
    temperature = GravityScalarRoughnessLength(FT)
    water_vapor = GravityScalarRoughnessLength(FT)
    return SimilarityScales(momentum, temperature, water_vapor)
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

function SimilarityTheoryTurbulentFluxes(FT::DataType = Float64;
                                         gravitational_acceleration = default_gravitational_acceleration,
                                         bulk_velocity_scale = nothing,
                                         von_karman_constant = convert(FT, 0.4),
                                         planetary_boundary_layer_height = convert(FT, 600),
                                         turbulent_prandtl_number = convert(FT, 1),
                                         similarity_functions = businger_similarity_functions(FT),
                                         thermodynamics_parameters = PATP(FT),
                                         water_vapor_saturation = ClasiusClapyeronSaturation(),
                                         water_mole_fraction = convert(FT, 0.98),
                                         roughness_lengths = default_roughness_lengths(FT),
                                         fields = nothing)

    return SimilarityTheoryTurbulentFluxes(convert(FT, gravitational_acceleration),
                                           convert(FT, von_karman_constant),
                                           convert(FT, planetary_boundary_layer_height),
                                           convert(FT, turbulent_prandtl_number),
                                           bulk_velocity_scale,
                                           similarity_functions,
                                           thermodynamics_parameters,
                                           water_vapor_saturation,
                                           water_mole_fraction,
                                           roughness_lengths,
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

@inline function seawater_saturation_specific_humidity(atmosphere_thermodynamics_parameters,
                                                       surface_temperature,
                                                       surface_salinity,
                                                       atmos_state,
                                                       water_mole_fraction,
                                                       water_vapor_saturation,
                                                       ::Liquid)

    ℂₐ = atmosphere_thermodynamics_parameters
    FT = eltype(ℂₐ)
    Tₛ = surface_temperature
    Sₛ = surface_salinity
    ρₛ = atmos_state.ρ # surface density -- should we extrapolate to obtain this?
    ρₛ = convert(FT, ρₛ)

    q★_H₂O = water_saturation_specific_humidity(water_vapor_saturation, ℂₐ, ρₛ, Tₛ)
    x_H₂O  = compute_water_mole_fraction(water_mole_fraction, Sₛ)

    # Return saturation specific humidity for salty seawater
    return q★_H₂O * x_H₂O
end

struct SalinityConstituent{FT}
    molar_mass :: FT
    mass_fraction :: FT
end

struct WaterMoleFraction{FT, C}
    water_molar_mass :: FT
    salinity_constituents :: C
end

function WaterMoleFraction(FT=Float64)
    water_molar_mass = convert(FT, 18.02)

    # TODO: find reference for these
    salinity_constituents = (
        chloride  = SalinityConstituent{FT}(35.45, 0.56),
        sodium    = SalinityConstituent{FT}(22.99, 0.31),
        sulfate   = SalinityConstituent{FT}(96.06, 0.08),
        magnesium = SalinityConstituent{FT}(24.31, 0.05),
    )

    return SeawaterComposition(water_molar_mass, salinity_constituents)
end

@inline compute_water_mole_fraction(x_H₂O::Number, S) = x_H₂O

@inline function compute_water_mole_fraction(wmf::WaterMoleFraction, S)
    # TODO: express the concept of "ocean_salinity_units"?
    s = S / 1000 # convert g/kg to concentration

    # Molecular weights
    μ_H₂O = wmf.water_molar_mass

    # Salinity constituents: Cl, Na, SO₄, Mg
    μ_Cl  = wmf.salinity_constituents.chloride.molar_mass
    μ_Na  = wmf.salinity_constituents.sodium.molar_mass
    μ_SO₄ = wmf.salinity_constituents.sulfate.molar_mass
    μ_Mg  = wmf.salinity_constituents.magnesium.molar_mass

    # Salinity constituent fractions
    ϵ_Cl  = wmf.salinity_constituents.chloride.mass_fraction
    ϵ_Na  = wmf.salinity_constituents.sodium.mass_fraction
    ϵ_SO₄ = wmf.salinity_constituents.sulfate.mass_fraction
    ϵ_Mg  = wmf.salinity_constituents.magnesium.mass_fraction

    α = μ_H₂O * (ϵ_Cl/μ_Cl + ϵ_Na/μ_Na  + ϵ_SO₄/μ_SO₄ + ϵ_Mg/μ_Mg)

    return (1 - s) / (1 - s + α * s)
end

#####
##### Struct that represents a 3-tuple of momentum, heat, and water vapor
#####

struct SimilarityScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

# Convenience default with water_vapor component = nothing
SimilarityScales(momentum, temperature) = SimilarityScales(momentum, temperature, nothing)

#####
##### Fixed-point iteration for roughness length
#####

@inline function compute_similarity_theory_fluxes(similarity_theory,
                                                  surface_state,
                                                  atmos_state,
                                                  thermodynamics_parameters,
                                                  gravitational_acceleration,
                                                  von_karman_constant)

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    Δh, Δu, Δv, Δθ, Δq = state_differences(ℂₐ, atmos_state, surface_state)
    differences = (; u=Δu, v=Δv, θ=Δθ, q=Δq, h=Δh)
    
    # Initial guess for the characteristic scales u★, θ★, q★.
    Σ★ = initial_guess(differences, 
                       similarity_theory,
                       gravitational_acceleration,
                       von_karman_constant, 
                       ℂₐ, atmos_state.ts)

    # The inital velocity scale assumes that
    # the gustiness velocity `uᴳ` is equal to 0.5 ms⁻¹. 
    # That will be refined later on.
    uτ = sqrt(Δu^2 + Δv^2 + 0.25)

    for _ in 1:10
        Σ★, uτ = refine_characteristic_scales(Σ★, uτ, 
                                              similarity_theory,
                                              surface_state,
                                              differences,
                                              thermodynamics_parameters,
                                              gravitational_acceleration,
                                              von_karman_constant)
    end

    u★ = Σ★.momentum
    θ★ = Σ★.temperature
    q★ = Σ★.water_vapor

    # u★² ≡ sqrt(τx² + τy²)
    τx = - u★^2 * Δu / sqrt(Δu^2 + Δv^2)
    τy = - u★^2 * Δv / sqrt(Δu^2 + Δv^2)

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

@inline function initial_guess(differences, 
                               similarity_theory,
                               gravitational_acceleration,
                               von_karman_constant, 
                               ℂₐ, 𝒬ₐ)

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q
    h  = differences.h

    g  = gravitational_acceleration
    ϰ  = von_karman_constant
    
    # Extract roughness lengths
    ℓu = similarity_theory.roughness_lengths.momentum
    zᵢ = similarity_theory.planetary_boundary_layer_height

    # assuming the initial gustiness is `0.5` ms⁻¹
    uᴳ = 0.5
    uτ = sqrt(Δu^2 + Δv^2 + uᴳ^2)

    # u10 at the reference ten meter height, assuming the initial roughness length is `1e-4` m
    u10 = uτ / log(h / 1e-4) * 11.5129 # log(10 / 1e-4) == 11.5129
    u★  = 0.035 * u10

    ℓu₀ = roughness_length(ℓu, u★, 𝒬ₐ, ℂₐ)

    # Initial neutral coefficients at 10 meter height
    Cuₙ  = (ϰ / log(10 / ℓu₀))^2
    Cθₙ  = 0.00115 / sqrt(Cuₙ)

    # Initial thermal roughness length
    ℓθ₀ = 10 / exp(ϰ / Cθₙ)

    # Neutral transfer coefficients at height `h`
    Cu = (ϰ / log(h / ℓu₀))^2
    Cθ = (ϰ / log(h / ℓθ₀))^2

    # Similarity functions from Businger et al. (1971)
    ψu = similarity_theory.similarity_functions.momentum
    ψθ = similarity_theory.similarity_functions.temperature
    ψq = similarity_theory.similarity_functions.water_vapor

    # Bulk Flux Richardson number
    Ribcu = - h / zᵢ / 0.004 / 1.2^3 # - h / zi / 0.004 / β^3
    b★    = buoyancy_scale(Δθ, Δq, 𝒬ₐ, ℂₐ, g)
    Ri    = - ifelse(b★ == 0, zero(b★), h / b★ / uτ^2)
    ζ     = ϰ * Cθ / Cu * Ribcu / (1 + Ri / Ribcu)
    
    L10   = ifelse(ζ == 0, zero(ζ), h / ζ)

    u★ = uτ * ϰ / (log(h / ℓu₀) - ψu(h / L10) + ψu(ℓu₀ / L10))
    θ★ = Δθ * ϰ / (log(h / ℓθ₀) - ψθ(h / L10) + ψθ(ℓθ₀ / L10))
    q★ = Δq * ϰ / (log(h / ℓθ₀) - ψq(h / L10) + ψq(ℓθ₀ / L10))
    
    return SimilarityScales(u★, θ★, q★)
end

struct Momentum end
struct Scalar end

struct SimilarityFunction{M, FT, C, D}
    a :: FT
    b :: FT
    c :: C
    d :: D

    SimilarityFunction{M}(a::FT, b::FT, c::C, d::D) where {M, FT, C, D} = new{M, FT, C, D}(a, b, c, d)
end

Adapt.adapt_structure(to, ψ::SimilarityFunction{M}) where M = SimilarityFunction{M}(ψ.a, ψ.b, ψ.c, ψ.d)

function businger_similarity_functions(FT = Float64)

    # Computed from Businger et al. (1971)
    # ψu = SimilarityFunction{Momentum}(4.7, 15.0, OneQuarter())
    # ψc = SimilarityFunction{Scalar}(6.35, 9.0, OneHalf())

    ψu = SimilarityFunction{Momentum}(0.7, 0.75, 5.0, 0.35)
    ψc = SimilarityFunction{Scalar}(6.35, 2/3, 0.35, 1.0)

    return SimilarityScales(ψu, ψc, ψc)
end

@inline function (ψ::SimilarityFunction{<:Momentum})(ζ)

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, 0.35 * ζ⁺)

    ψ_stable = - 0.7 * ζ⁺ - 3 / 4 * (ζ⁺ - 5 / 0.35) * exp(-dζ) - 3 / 4 * 5 / 0.35
    
    fₘ = sqrt(sqrt(1 - 15 * ζ⁻))
    ψ_unstable_1 = log((1 + fₘ)^2 * (1 + fₘ^2) / 8) - 2 * atan(fₘ) + π / 2;

    fₘ = cbrt(1 - 10.15 * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₘ + fₘ^2) / 3) - sqrt(3) * atan((1 + 2fₘ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

@inline function (ψ::SimilarityFunction{<:Scalar})(ζ)

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, 0.35 * ζ⁺)

    ψ_stable = - (4 * ζ⁺ / 3)^(3 / 2) - 2 / 3 * (ζ⁺ - 14.28) * exp(-dζ) - 8.525
    
    fₕ = sqrt(1 - 15 * ζ⁻)
    ψ_unstable_1 = 2 * log((1 + fₕ) / 2) 

    fₕ = cbrt(1 - 34.15 * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₕ + fₕ^2) / 3) - sqrt(3) * atan((1 + 2fₕ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

struct OneQuarter end
struct OneHalf end

import Base: ^
@inline ^(x, ::OneQuarter) = sqrt(sqrt(x))
@inline ^(x, ::OneHalf) = sqrt(x)

# The M-O characteristic length is calculated as
#  L★ = - u★² / (κ ⋅ b★)
# where b★ is the characteristic buoyancy scale calculated from this function
@inline function buoyancy_scale(θ★, q★, 𝒬, ℂ, g)
    𝒯₀ = AtmosphericThermodynamics.virtual_temperature(ℂ, 𝒬)
    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬)

    ε = AtmosphericThermodynamics.Parameters.molmass_ratio(ℂ)
    δ = ε - 1 # typically equal to 0.608

    # Where does this come from? Probably Fairell et al. 1996, 
    b★ = g / 𝒯₀ * (θ★ * (1 + δ * q₀) + δ * θ₀ * q★)

    return b★
end

@inline function state_differences(ℂ, 𝒰₁, 𝒰₀)
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
    Δθ = θ₁ - θ₀

    q₁ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₁)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₀)
    Δq = q₁ - q₀

    return Δh, Δu, Δv, Δθ, Δq
end

@inline roughness_length(ℓ, u★, args...)     = ℓ(u★, args...)
@inline roughness_length(ℓ::Number, args...) = ℓ

@inline function refine_characteristic_scales(estimated_characteristic_scales, 
                                              velocity_scale,
                                              similarity_theory,
                                              surface_state,
                                              differences,
                                              thermodynamics_parameters,
                                              gravitational_acceleration,
                                              von_karman_constant)

    # "initial" scales because we will recompute them
    u★ = estimated_characteristic_scales.momentum
    θ★ = estimated_characteristic_scales.temperature
    q★ = estimated_characteristic_scales.water_vapor
    uτ = velocity_scale

    # Similarity functions from Businger et al. (1971)
    ψu = similarity_theory.similarity_functions.momentum
    ψθ = similarity_theory.similarity_functions.temperature
    ψq = similarity_theory.similarity_functions.water_vapor

    # Extract roughness lengths
    ℓu = similarity_theory.roughness_lengths.momentum
    ℓθ = similarity_theory.roughness_lengths.temperature
    ℓq = similarity_theory.roughness_lengths.water_vapor

    h  = differences.h
    ϰ  = von_karman_constant
    ℂ  = thermodynamics_parameters
    g  = gravitational_acceleration
    𝒬ₒ = surface_state.ts # thermodynamic state
    zᵢ = similarity_theory.planetary_boundary_layer_height

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    b★ = buoyancy_scale(θ★, q★, 𝒬ₒ, ℂ, g)

    # Monin-Obhukov characteristic length scale and non-dimensional height
    L★ = ifelse(b★ == 0, zero(b★), - u★^2 / (ϰ * b★))

    # Compute roughness length scales
    ℓu₀ = roughness_length(ℓu, u★, 𝒬ₒ, ℂ)
    ℓq₀ = roughness_length(ℓq, ℓu₀, u★, 𝒬ₒ, ℂ)
    ℓθ₀ = roughness_length(ℓθ, ℓu₀, u★, 𝒬ₒ, ℂ)

    # Transfer coefficients at height `h`
    Cu = ϰ / (log(h / ℓu₀) - ψu(h / L★) + ψu(ℓu₀ / L★))
    Cθ = ϰ / (log(h / ℓq₀) - ψθ(h / L★) + ψθ(ℓq₀ / L★))
    Cq = ϰ / (log(h / ℓθ₀) - ψq(h / L★) + ψq(ℓθ₀ / L★))

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    # u★ including gustiness
    u★ = Cu * uτ
    θ★ = Cθ * Δθ
    q★ = Cq * Δq

    # Dissipation characteristic scale for gustiness
    ε★ = - u★ * b★
    uᴳ = 1.2 * cbrt(ε★ * zᵢ)

    # New velocity difference accounting for gustiness
    uτ = sqrt(Δu^2 + Δv^2 + uᴳ^2)

    return SimilarityScales(u★, θ★, q★), uτ
end

struct GravityMomentumRoughnessLength{FT, V}
    gravitational_acceleration :: FT
    air_kinematic_viscosity :: V
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
    maximum_roughness_length :: FT
end

struct GravityScalarRoughnessLength{FT, V, R}
    air_kinematic_viscosity :: V
    reynolds_number_scaling_function :: R
    maximum_roughness_length :: FT
end

# Empirical fit of the scalar roughness length with roughness Reynolds number `R★ = u★ / ν`
# Edson et al. (2013), equation (28)
@inline empirical_scaling_function(R★ :: FT, args...) where FT = 
        ifelse(R★ == 0, FT(0), convert(FT, 5.85e-5 / R★ ^ 0.72))

# Assumes that θ comes in in Kelvin
@inline function temperature_dependent_viscosity(θ :: FT) where FT 
    T = convert(FT, θ - celsius_to_kelvin)
    ν = convert(FT, 1.326e-5 * (1 + 6.542e-3 * T + 8.301e-6 * T^2 - 4.84e-9 * T^3))
    
    return ν
end

function GravityScalarRoughnessLength(FT=Float64;
                                      air_kinematic_viscosity = temperature_dependent_viscosity,
                                      reynolds_number_scaling_function = empirical_scaling_function,
                                      maximum_roughness_length = 1.6e-4) # Values from COARE3.6

    return GravityScalarRoughnessLength(air_kinematic_viscosity,
                                        reynolds_number_scaling_function,
                                        convert(FT, maximum_roughness_length))
end

function GravityMomentumRoughnessLength(FT=Float64;
                                        gravitational_acceleration = default_gravitational_acceleration,
                                        maximum_roughness_length = Inf, # An estimate?
                                        air_kinematic_viscosity = temperature_dependent_viscosity,
                                        gravity_wave_parameter = 0.011,
                                        laminar_parameter = 0.11)

    return GravityMomentumRoughnessLength(convert(FT, gravitational_acceleration),
                                          air_kinematic_viscosity,
                                          convert(FT, gravity_wave_parameter),
                                          convert(FT, laminar_parameter),
                                          convert(FT, maximum_roughness_length))
end

# Momentum roughness length should be different from scalar roughness length.
# Apparently temperature and water vapor can be considered the same (Edison et al 2013)
@inline function roughness_length(ℓ::GravityMomentumRoughnessLength{FT}, u★, 𝒬, ℂ) where FT
    g  = ℓ.gravitational_acceleration
    α  = ℓ.gravity_wave_parameter
    β  = ℓ.laminar_parameter
    ℓm = ℓ.maximum_roughness_length

    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬)
    ν  = ℓ.air_kinematic_viscosity(θ₀)

    # We need to prevent `Inf` that pops up when `u★ == 0`.
    # For this reason, if `u★ == 0` we prescribe the roughness length to be
    # equal to a `maximum` roughness length
    ℓᴿ = ifelse(u★ == 0, ℓm, β * ν / u★) 
    
    return min(α * u★^2 / g + ℓᴿ, ℓm)
end

# This, for example is what is implemented in COARE 3.6
@inline function roughness_length(ℓ::GravityScalarRoughnessLength{FT}, ℓu, u★, 𝒬, ℂ) where FT
    ℓm = ℓ.maximum_roughness_length
    
    scaling_function = ℓ.reynolds_number_scaling_function

    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬)
    ν  = ℓ.air_kinematic_viscosity(θ₀)

    # Roughness Reynolds number
    R★ = ℓu * u★ / ν

    # implementation of scalar roughness length
    ℓq = scaling_function(R★, ℓu, u★, ν)

    return min(ℓq, ℓm) 
end