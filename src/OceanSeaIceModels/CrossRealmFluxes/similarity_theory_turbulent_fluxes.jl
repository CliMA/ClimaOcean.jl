using Oceananigans.Utils: prettysummary
using Oceananigans.Grids: AbstractGrid

using Thermodynamics: Liquid
using SurfaceFluxes.Parameters: SurfaceFluxesParameters, AbstractSurfaceFluxesParameters
using SurfaceFluxes.UniversalFunctions: BusingerParams

using ..PrescribedAtmospheres: PrescribedAtmosphereThermodynamicsParameters

import Thermodynamics as AtmosphericThermodynamics

import SurfaceFluxes.Parameters:
    thermodynamics_params,
    uf_params,
    von_karman_const,
    universal_func_type,
    grav

#####
##### Bulk turbulent fluxes based on similarity theory
#####

struct SimilarityTheoryTurbulentFluxes{FT, ΔU, UF, TP, S, W, F} <: AbstractSurfaceFluxesParameters
    gravitational_acceleration :: FT
    von_karman_constant :: FT
    bulk_velocity_scale :: ΔU
    universal_function :: UF
    thermodynamics_parameters :: TP
    water_vapor_saturation :: S
    water_mole_fraction :: W
    fields :: F
end

const STTF = SimilarityTheoryTurbulentFluxes
thermodynamics_params(fluxes::STTF) = fluxes.thermodynamics_parameters
universal_func_type(fluxes::STTF)   = universal_func_type(typeof(fluxes.universal_function))
uf_params(fluxes::STTF)             = fluxes.universal_function
von_karman_const(fluxes::STTF)      = fluxes.von_karman_constant
grav(fluxes::STTF)                  = fluxes.gravitational_acceleration

Base.summary(::SimilarityTheoryTurbulentFluxes{FT}) where FT = "SimilarityTheoryTurbulentFluxes{$FT}"

struct ClasiusClapyeronSaturation end
 
@inline function water_saturation_specific_humidity(::ClasiusClapyeronSaturation, ℂₐ, ρₛ, Tₛ)
    p★ = AtmosphericThermodynamics.saturation_vapor_pressure(ℂₐ, Tₛ, Liquid())
    q★ = AtmosphericThermodynamics.q_vap_saturation_from_density(ℂₐ, Tₛ, ρₛ, p★)
    return q★
end

Base.@kwdef struct LargeYeagerSaturation{FT}
    c₁ :: FT = 640380.0 # kg m⁻³
    c₂ :: FT = 5107.4   # K
end

const LYS = LargeYeagerSaturation
@inline water_saturation_specific_humidity(lys::LYS, ℂₐ, ρₛ, Tₛ) = lys.c₁ * exp(-lys.c₂ / Tₛ) / ρₛ

function Base.show(io::IO, fluxes::SimilarityTheoryTurbulentFluxes)
    print(io, summary(fluxes), '\n',
          "├── gravitational_acceleration: ",   prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",          prettysummary(fluxes.von_karman_constant), '\n',
          "├── bulk_velocity_scale: ",          summary(fluxes.bulk_velocity_scale), '\n',
          "├── universal_function: ",           summary(fluxes.universal_function), '\n',
          "├── water_mole_fraction: ",          summary(fluxes.water_mole_fraction), '\n',
          "├── water_vapor_saturation: ",       summary(fluxes.water_vapor_saturation), '\n',
          "└── thermodynamics_parameters: ",    summary(fluxes.thermodynamics_parameters))
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

function SimilarityTheoryTurbulentFluxes(FT::DataType = Float64;
                                         gravitational_acceleration = convert(FT, 9.80665),
                                         bulk_velocity_scale = nothing,
                                         von_karman_constant = convert(FT, 0.4),
                                         universal_function = default_universal_function_parameters(FT),
                                         thermodynamics_parameters = PATP(FT),
                                         water_vapor_saturation = ClasiusClapyeronSaturation(),
                                         water_mole_fraction = convert(FT, 0.98),
                                         fields = nothing)

    return SimilarityTheoryTurbulentFluxes(gravitational_acceleration,
                                           von_karman_constant,
                                           bulk_velocity_scale,
                                           universal_function,
                                           thermodynamics_parameters,
                                           water_vapor_saturation,
                                           water_mole_fraction,
                                           fields)
end

function SimilarityTheoryTurbulentFluxes(grid::AbstractGrid; kw...)
    evaporation = Field{Center, Center, Nothing}(grid)
    latent_heat_flux = Field{Center, Center, Nothing}(grid)
    sensible_heat_flux = Field{Center, Center, Nothing}(grid)

    fields = (; latent_heat_flux, sensible_heat_flux, evaporation)

    return SimilarityTheoryTurbulentFluxes(eltype(grid); kw..., fields)
end

@inline update_turbulent_flux_fields!(::Nothing, args...) = nothing

@inline clip(x::FT) = max(zero(FT), x)

@inline function update_turbulent_flux_fields!(fields, i, j, grid, conditions, ρᶠ)
    Qᵉ = fields.latent_heat_flux
    Qᶜ = fields.sensible_heat_flux
    E = fields.evaporation
    @inbounds begin
        # +0: cooling, -0: heating
        Qᵉ[i, j, 1] = clip(conditions.lhf)
        Qᶜ[i, j, 1] = conditions.shf

        # "Salt flux" has the opposite sign of "freshwater flux"
        Eᵢ = - conditions.evaporation / ρᶠ # convert to volume flux
        E[i, j, 1] = clip(Eᵢ)
    end
    return nothing
end

# See SurfaceFluxes.jl for other parameter set options.
default_universal_function_parameters(FT=Float64) = BusingerParams{FT}(Pr_0 = convert(FT, 0.74),
                                                                       a_m  = convert(FT, 4.7),
                                                                       a_h  = convert(FT, 4.7),
                                                                       ζ_a  = convert(FT, 2.5),
                                                                       γ    = convert(FT, 4.42))

function seawater_saturation_specific_humidity(atmosphere_thermodynamics_parameters,
                                               surface_temperature,
                                               surface_salinity,
                                               atmos_state,
                                               water_mole_fraction,
                                               water_vapor_saturation,
                                               ::Liquid)

    ℂₐ = atmosphere_thermodynamics_parameters
    Tₛ = surface_temperature
    Sₛ = surface_salinity
    ρₛ = atmos_state.ρ # surface density -- should we extrapolate to obtain this?

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

