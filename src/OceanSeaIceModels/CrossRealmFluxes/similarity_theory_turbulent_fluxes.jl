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
          "├── gravitational_acceleration: ",   prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",          prettysummary(fluxes.von_karman_constant), '\n',
          "├── bulk_velocity_scale: ",          summary(fluxes.bulk_velocity_scale), '\n',
          "├── similarity_function: ",          summary(fluxes.similarity_function), '\n',
          "├── water_mole_fraction: ",          summary(fluxes.water_mole_fraction), '\n',
          "├── water_vapor_saturation: ",       summary(fluxes.water_vapor_saturation), '\n',
          "└── thermodynamics_parameters: ",    summary(fluxes.thermodynamics_parameters))
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
                                         similarity_functions = businger_similarity_functions(FT),
                                         thermodynamics_parameters = PATP(FT),
                                         water_vapor_saturation = ClasiusClapyeronSaturation(),
                                         water_mole_fraction = convert(FT, 0.98),
                                         roughness_lengths = default_roughness_lengths(FT),
                                         fields = nothing)

    return SimilarityTheoryTurbulentFluxes(convert(FT, gravitational_acceleration),
                                           convert(FT, von_karman_constant),
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

@inline function compute_similarity_theory_fluxes(roughness_lengths,
                                                  similarity_functions,
                                                  surface_state,
                                                  atmos_state,
                                                  thermodynamics_parameters,
                                                  gravitational_acceleration,
                                                  von_karman_constant,
                                                  Σ₀ = SimilarityScales(1e-3, 1e-3, 1e-3))

    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    Δh, Δu, Δv, Δθ, Δq = state_differences(ℂₐ, atmos_state, surface_state)
    differences = (; u=Δu, v=Δv, θ=Δθ, q=Δq, h=Δh)

    # Solve for the characteristic scales u★, θ★, q★, and thus for fluxes.
    Σ★ = Σ₀

    for _ in 1:10
        Σ★ = refine_characteristic_scales(Σ★, 
                                          roughness_lengths,
                                          similarity_functions, 
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

struct Momentum end
struct Scalar end

struct SimilarityFunction{M, FT, C}
    a :: FT
    b :: FT
    c :: C

    SimilarityFunction{M}(a::FT, b::FT, c::C) where {M, FT, C} = new{M, FT, C}(a, b, c)
end

Adapt.adapt_structure(to, ψ::SimilarityFunction{M}) where M = SimilarityFunction{M}(ψ.a, ψ.b, ψ.c)

function businger_similarity_functions(FT = Float64)

    # Computed from Businger et al. (1971)
    ψu = SimilarityFunction{Momentum}(4.7, 15.0, OneQuarter())
    ψc = SimilarityFunction{Scalar}(6.35, 9.0, OneHalf())

    return SimilarityScales(ψu, ψc, ψc)
end

# This seems to come from "SURFACE FLUXES FOR PRACTITIONERS OF GLOBAL OCEAN DATA ASSIMILATION"
# Of William Large, but a couple of coefficients and signs are off.
# Also in that paper momentum and scalar stability functions are different, here they are the same??
# Fairell et al implement a different formulation with a "convective" and "stable" stability function
@inline function (ψ::SimilarityFunction{<:Momentum})(ζ)
    a = ψ.a
    b = ψ.b
    c = ψ.c

    ζ⁻ = min(zero(ζ), ζ)
    fₘ = (1 - b * ζ⁻)^c

    ψ_unstable = log((1 + fₘ)^2 * (1 + fₘ^2) / 8) - 2 * atan(fₘ) + π / 2
    ψ_stable   = - a * ζ

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

@inline function (ψ::SimilarityFunction{<:Scalar})(ζ)
    a = ψ.a
    b = ψ.b
    c = ψ.c

    ζ⁻ = min(zero(ζ), ζ)
    fₕ = (1 - b * ζ⁻)^c

    ψ_unstable = 2 * log((1 + fₕ^2) / 2) 
    ψ_stable   = - a * ζ

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

struct OneQuarter end
struct OneHalf end

import Base: ^
@inline ^(x, ::OneQuarter) = sqrt(sqrt(x))
@inline ^(x, ::OneHalf) = sqrt(x)

@inline function bulk_factor(ψ, h, ℓ, L★)

    # Non-dimensional height in Obukhov length units
    ζ  = ifelse(L★ == 0, zero(h), h / L★) 

    # Non-dimensional roughness height in Obukhov length units
    ζᵣ = ifelse(L★ == 0, zero(h), ℓ / L★) 

    χ⁻¹ = log(h / ℓ) - ψ(ζ) + ψ(ζᵣ)
    
    return ifelse(χ⁻¹ == 0, zero(h), 1 / χ⁻¹)
end

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

@inline roughness_length(ℓ, Σ★) = ℓ(Σ★)
@inline roughness_length(ℓ::Number, Σ★) = ℓ

@inline roughness_length(ℓ, ℓu, Σ★) = ℓ(Σ★)
@inline roughness_length(ℓ::Number, ℓu, Σ★) = ℓ

@inline function refine_characteristic_scales(estimated_characteristic_scales,
                                              roughness_lengths,
                                              similarity_functions,
                                              surface_state,
                                              differences,
                                              thermodynamics_parameters,
                                              gravitational_acceleration,
                                              von_karman_constant)

    # "initial" scales because we will recompute them
    u★ = estimated_characteristic_scales.momentum
    θ★ = estimated_characteristic_scales.temperature
    q★ = estimated_characteristic_scales.water_vapor
    Σ★ = estimated_characteristic_scales

    # Similarity functions from Businger et al. (1971)
    ψu = similarity_functions.momentum
    ψθ = similarity_functions.temperature
    ψq = similarity_functions.water_vapor

    # Extract roughness lengths
    ℓu = roughness_lengths.momentum
    ℓθ = roughness_lengths.temperature
    ℓq = roughness_lengths.water_vapor

    h = differences.h
    ϰ = von_karman_constant
    
    # Compute roughness length scales
    ℓu₀ = roughness_length(ℓu, Σ★)
    ℓq₀ = roughness_length(ℓq, ℓu₀, Σ★)
    ℓθ₀ = roughness_length(ℓθ, ℓu₀, Σ★)

    # Compute Monin-Obukhov length scale depending on a `buoyancy flux`
    ℂ = thermodynamics_parameters
    g = gravitational_acceleration
    𝒬ₒ = surface_state.ts # thermodynamic state
    b★ = buoyancy_scale(θ★, q★, 𝒬ₒ, ℂ, g)

    # Monin-Obhukov characteristic length scale
    L★ = ifelse(b★ == 0, zero(b★), - u★^2 / (ϰ * b★))

    χu = bulk_factor(ψu, h, ℓu₀, L★)
    χθ = bulk_factor(ψθ, h, ℓθ₀, L★)
    χq = bulk_factor(ψq, h, ℓq₀, L★)

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    # Maybe we should add gustiness here?
    u★ = ϰ * χu * sqrt(Δu^2 + Δv^2) 
    θ★ = ϰ * χθ * Δθ
    q★ = ϰ * χq * Δq

    return SimilarityScales(u★, θ★, q★)
end

struct GravityMomentumRoughnessLength{FT}
    gravitational_acceleration :: FT
    air_kinematic_viscosity :: FT
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
    maximum_roughness_length :: FT
end

struct GravityScalarRoughnessLength{FT, R}
    air_kinematic_viscosity :: FT
    reynolds_number_scaling_function :: R
    maximum_roughness_length :: FT
end

# Empirical fit of the scalar roughness length with roughness Reynolds number `R★ = u★ / ν`
# Edson et al. (2013), equation (28)
@inline empirical_scaling_function(R★ :: FT, args...) where FT = 
        ifelse(R★ == 0, FT(0), convert(FT, 5.85e-5 / R★ ^ 0.76))

# Brusser - Garrat scaling of the scalar roughness length with roughness number 
# Edson et al. (2013), equation (29)
@inline brusser_garrat_scaling_function(R★ :: FT, ℓu, args...) where FT = 
        convert(FT, ℓu * exp(2 - 2.28 * sqrt(sqrt(R★))))


function GravityScalarRoughnessLength(FT=Float64;
                                      air_kinematic_viscosity = 1.5e-5,
                                      reynolds_number_scaling_function = brusser_garrat_scaling_function,
                                      maximum_roughness_length = 1.6e-4) # Values from COARE3.6

    return GravityScalarRoughnessLength(convert(FT, air_kinematic_viscosity),
                                        reynolds_number_scaling_function,
                                        convert(FT, maximum_roughness_length))
end

function GravityMomentumRoughnessLength(FT=Float64;
                                    gravitational_acceleration = default_gravitational_acceleration,
                                    maximum_roughness_length = 5e-3, # An estimate?
                                    air_kinematic_viscosity = 1.5e-5,
                                    gravity_wave_parameter = 0.011,
                                    laminar_parameter = 0.11)

    return GravityMomentumRoughnessLength(convert(FT, gravitational_acceleration),
                                          convert(FT, air_kinematic_viscosity),
                                          convert(FT, gravity_wave_parameter),
                                          convert(FT, laminar_parameter),
                                          convert(FT, maximum_roughness_length))
end

# Momentum roughness length should be different from scalar roughness length.
# Apparently temperature and water vapor can be considered the same (Edison et al 2013)
@inline function roughness_length(ℓ::GravityMomentumRoughnessLength{FT}, Σ★) where FT
    u★ = Σ★.momentum
    g  = ℓ.gravitational_acceleration
    ν  = ℓ.air_kinematic_viscosity
    α  = ℓ.gravity_wave_parameter
    β  = ℓ.laminar_parameter
    ℓm = ℓ.maximum_roughness_length

    # We need to prevent `Inf` that pops up when `u★ == 0`.
    # For this reason, if `u★ == 0` we prescribe the roughness length to be
    # equal to a `maximum` roughness length
    ℓᴿ = ifelse(u★ == 0, ℓm, β * ν / u★) 
    
    return min(α * u★^2 / g + ℓᴿ, ℓm)
end

# This, for example is what is implemented in COARE 3.6
@inline function roughness_length(ℓ::GravityScalarRoughnessLength{FT}, ℓu, Σ★) where FT
    u★ = Σ★.momentum
    ν  = ℓ.air_kinematic_viscosity
    ℓm = ℓ.maximum_roughness_length
    
    scaling_function = ℓ.reynolds_number_scaling_function

    # Roughness Reynolds number
    R★ = ℓu * u★ / ν

    # implementation of scalar roughness length
    ℓq = scaling_function(R★, ℓu, u★, ν)

    return min(ℓq, ℓm) 
end