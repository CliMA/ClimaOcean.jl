using Oceananigans.Utils: prettysummary
using Oceananigans.Grids: AbstractGrid

using Adapt
using Thermodynamics: Liquid
using SurfaceFluxes.Parameters: SurfaceFluxesParameters, AbstractSurfaceFluxesParameters
using SurfaceFluxes.UniversalFunctions: BusingerParams, BusingerType

using ..PrescribedAtmospheres: PrescribedAtmosphereThermodynamicsParameters

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
                                                                          adapt(to, fluxes.bulk_velocity_scale),
                                                                          adapt(to, fluxes.similarity_function),
                                                                          adapt(to, fluxes.thermodynamics_parameters),
                                                                          adapt(to, fluxes.water_vapor_saturation),
                                                                          adapt(to, fluxes.water_mole_fraction),
                                                                          adapt(to, fluxes.roughness_lengths),
                                                                          adapt(to, fluxes.fields))

Base.summary(::SimilarityTheoryTurbulentFluxes{FT}) where FT = "SimilarityTheoryTurbulentFluxes{$FT}"

struct ClasiusClapyeronSaturation end
 
@inline function water_saturation_specific_humidity(::ClasiusClapyeronSaturation, ℂₐ, ρₛ, Tₛ)
    p★ = AtmosphericThermodynamics.saturation_vapor_pressure(ℂₐ, Tₛ, Liquid())
    q★ = AtmosphericThermodynamics.q_vap_saturation_from_density(ℂₐ, Tₛ, ρₛ, p★)
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
    momentum = convert(FT, 1e-4)
    heat     = convert(FT, 1e-4)
    return SimilarityScales(momentum, heat)
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
                                         # roughness_lengths = default_roughness_lengths(FT),
                                         roughness_lengths = SimilarityScales(1e-3, 1e-3, 1e-3),
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

# See SurfaceFluxes.jl for other parameter set options.
default_businger_parameters(FT=Float64) = BusingerParams{FT}(Pr_0 = convert(FT, 0.74),
                                                             a_m  = convert(FT, 4.7),
                                                             a_h  = convert(FT, 4.7),
                                                             ζ_a  = convert(FT, 2.5),
                                                             γ    = convert(FT, 4.42))

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

#=
struct SimilarityFunction{FT}
    a :: FT
    b :: FT
    c :: FT
end

struct GravityWaveRoughnessLength{FT}
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
    air_kinematic_viscosity :: FT
end

struct AtmosphericState{Q, T, U, V}
    q :: Q
    θ :: T
    u :: U
    v :: V
end

AtmosphericState(q, θ, u) = AtmosphericState(q, θ, u, nothing)

@inline function (ψ::SimilarityFunction)(Ri)
    a = ψ.a
    b = ψ.b
    c = ψ.c

    ϕ⁻¹ = (1 - b * Ri)^c
    ψ_unstable = log((1 + ϕ⁻¹)^2 * (1 + ϕ⁻¹^2) / 8) - 2 * atan(ϕ⁻¹) + π/2
    ψ_stable = - a * Ri
    return ifelse(Ri < 0, ψ_unstable, ψ_stable)
end

@inline similarity_scale(ψ, h, ℓ, Ri) = 1 / (log(h/ℓ) - ψ(Ri) + ψ(ℓ * Ri / h))

function buoyancy_scale(θ★, q★, surface_state, parameters)
    θ★ = fluxes.θ
    q★ = fluxes.q
    𝒯₀ = virtual_temperature(parameters, surface_state)
    q₀ = surface_state.q
    θ₀ = surface_state.θ
    r = parameters.molar_mass_ratio
    g = parameters.gravitational_acceleration
    δ = r - 1
    b★ = g / 𝒯₀ * (θ★ * (1 + δ * q₀) + δ * θ₀ * q★)
    return b★
end

function fixed_point_fluxes(u★, θ★, q★,
                            surface_state,
                            inner_length_scales,
                            universal_function,
                            parameters)

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    ϰ = parameters.von_karman_constant
    f = universal_function

    b★ = buoyancy_scale(θ★, q★, surface_state, parameters)
    Riₕ = - ϰ * h * b★ / u★^2

    ℓu = inner_length_scales.u(u★)
    ℓθ = inner_length_scales.θ(u★)
    ℓq = inner_length_scales.q(u★)

    χu = momentum_flux_scale(f, h, ℓu, Riₕ)
    χθ =   tracer_flux_scale(f, h, ℓθ, Riₕ)
    χq =   tracer_flux_scale(f, h, ℓq, Riₕ)

    u★ = ϰ * χu * sqrt(Δu^2 + Δv^2)
    θ★ = ϰ * χθ * Δθ
    q★ = ϰ * χq * Δq

    return u★, θ★, q★
end

function GravityWaveRoughnessLengths(FT=Float64;
                                     gravity_wave_parameter = 0.011,
                                     laminar_parameter = 0.11,
                                     air_kinematic_viscosity=1.5e-5)

    return GravityWaveRoughnessLengths(convert(FT, gravity_wave_parameter),
                                       convert(FT, laminar_parameter),
                                       convert(FT, air_kinematic_viscosity))
end

@inline function compute_turbulent_surface_fluxes(similarity_function::BusingerParams,
                                                  roughness_lengths,
                                                  atmos_state,
                                                  ocean_state)

    ℓu = roughness_lengths.momentum
    ℓθ = roughness_lengths.heat
    ℓq = roughness_lengths.water_vapor
                                                    

    fluxes = (;
        latent_heat_flux         = conditions.lhf,
        sensible_heat_flux       = conditions.shf,
        freshwater_flux          = conditions.evaporation,
        zonal_momentum_flux      = conditions.ρτxz,
        meridional_momentum_flux = conditions.ρτyz,
    )

@inline function compute_turbulent_surface_fluxes(similarity_function::BusingerParams,
                                                  roughness_lengths::SimplifiedRoughnessLengths,
                                                  atmos_state,
                                                  ocean_state)

    # Solve for the surface fluxes with initial roughness length guess
    Uᵍ = zero(grid) # gustiness
    β = one(grid)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(atmos_state, ocean_State,
                                      roughness_lengths.momentum,
                                      roughness_lengths.heat
                                      Uᵍ, β)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)

    fluxes = (;
        latent_heat_flux         = conditions.lhf,
        sensible_heat_flux       = conditions.shf,
        freshwater_flux          = conditions.evaporation,
        zonal_momentum_flux      = conditions.ρτxz,
        meridional_momentum_flux = conditions.ρτyz,
    )

    return fluxes
end


@inline function compute_turbulent_surface_fluxes(roughness_lengths::GravityWaveRoughnessLengths,
                                                  atmos_state,
                                                  ocean_state)

    # Solve for the surface fluxes with initial roughness length guess
    Uᵍ = zero(grid) # gustiness
    β = one(grid)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(atmos_state, ocean_State,
                                      roughness_lengths.momentum,
                                      roughness_lengths.heat
                                      Uᵍ, β)

    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)

    fluxes = (;
        latent_heat_flux         = conditions.lhf,
        sensible_heat_flux       = conditions.shf,
        freshwater_flux          = conditions.evaporation,
        zonal_momentum_flux      = conditions.ρτxz,
        meridional_momentum_flux = conditions.ρτyz,
    )

    return fluxes
end

=#

