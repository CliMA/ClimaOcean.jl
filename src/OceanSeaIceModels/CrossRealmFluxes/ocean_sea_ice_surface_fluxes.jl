using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
using ClimaSeaIce.SlabSeaIceModels: SlabSeaIceModel

#####
##### Container for organizing information related to fluxes
#####

struct OceanSeaIceSurfaceFluxes{C, R, T, P}
    total :: C
    emitted_radiation :: R
    turbulent :: T
    prescribed :: P
end

Base.summary(crf::OceanSeaIceSurfaceFluxes) = "OceanSeaIceSurfaceFluxes"
Base.show(io::IO, crf::OceanSeaIceSurfaceFluxes) = print(io, summary(crf))

function OceanSeaIceSurfaceFluxes(ocean, sea_ice=nothing;
                                  atmosphere = nothing,
                                  surface_radiation = nothing)

    FT = eltype(ocean.model.grid)
    turbulent_fluxes = SimilarityTheoryTurbulentFluxes(FT)
    prescribed_fluxes = nothing

    return OceanSeaIceSurfaceFluxes(nothing,
                                    surface_radiation,
                                    turbulent_fluxes,
                                    prescribed_fluxes)
end

#=
function default_atmosphere_ocean_fluxes(FT=Float64, tracers=tuple(:S))
    # Note: we are constantly coping with the fact that the ocean is ᵒC.
    ocean_reference_temperature = 273.15
    momentum_transfer_coefficient = 5e-3
    evaporation_transfer_coefficient = 1e-3
    sensible_heat_transfer_coefficient = 2e-3
    vaporization_enthalpy  = 2.5e-3

    momentum_transfer_coefficient      = convert(FT, momentum_transfer_coefficient)
    evaporation_transfer_coefficient   = convert(FT, evaporation_transfer_coefficient)
    sensible_heat_transfer_coefficient = convert(FT, sensible_heat_transfer_coefficient)
    vaporization_enthalpy              = convert(FT, vaporization_enthalpy)

    τˣ = BulkFormula(RelativeUVelocity(), momentum_transfer_coefficient)
    τʸ = BulkFormula(RelativeVVelocity(), momentum_transfer_coefficient)
    momentum_flux_formulae = (u=τˣ, v=τʸ)

    # Note: reference temperature comes in here
    water_specific_humidity_difference = SpecificHumidity(FT)
    evaporation = nothing #BulkFormula(SpecificHumidity(FT), evaporation_transfer_coefficient)
    tracer_flux_formulae = (; S = evaporation)

    latent_heat_difference = LatentHeat(specific_humidity_difference = water_specific_humidity_difference; vaporization_enthalpy)
    latent_heat_formula    = nothing #BulkFormula(latent_heat_difference,  evaporation_transfer_coefficient)

    sensible_heat_difference = SensibleHeat(FT; ocean_reference_temperature)
    sensible_heat_formula = BulkFormula(sensible_heat_difference, sensible_heat_transfer_coefficient)

    heat_flux_formulae = (sensible_heat_formula, latent_heat_formula)

    return CrossRealmFluxes(momentum = momentum_flux_formulae,
                            heat = heat_flux_formulae,
                            tracers = tracer_flux_formulae)
end
=#

#=
#####
##### Bulk formula
#####

"""
    BulkFormula(air_sea_difference, transfer_coefficient)

The basic structure of a flux `J` computed by a bulk formula is:

```math
J = - ρₐ * C * Δc * ΔU
```

where `ρₐ` is the density of air, `C` is the `transfer_coefficient`,
`Δc` is the air_sea_difference, and `ΔU` is the bulk velocity scale.
"""
struct BulkFormula{F, CD}
    air_sea_difference :: F
    transfer_coefficient :: CD
end

@inline function cross_realm_flux(i, j, grid, time, formula::BulkFormula, ΔU, atmos_state, ocean_state)
    ρₐ = stateindex(atmos_state.ρ, i, j, 1, time)
    C = formula.transfer_coefficient
    Δc = air_sea_difference(i, j, grid, time, formula.air_sea_difference, atmos_state, ocean_state)

    # Note the sign convention, which corresponds to positive upward fluxes:
    return - ρₐ * C * Δc * ΔU
end

@inline cross_realm_flux(i, j, grid, time, ::Nothing,        args...) = zero(grid)
@inline cross_realm_flux(i, j, grid, time, a::AbstractArray, args...) = stateindex(a, i, j, 1, time)
@inline cross_realm_flux(i, j, grid, time, nt::NamedTuple,   args...) = cross_realm_flux(i, j, grid, time, values(nt), args...)

@inline cross_realm_flux(i, j, grid, time, flux_tuple::Tuple{<:Any, <:Any}, args...) =
    cross_realm_flux(i, j, grid, time, flux_tuple[1], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[2], args...)

@inline cross_realm_flux(i, j, grid, time, flux_tuple::Tuple{<:Any, <:Any, <:Any}, args...) =
    cross_realm_flux(i, j, grid, time, flux_tuple[1], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[2], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[3], args...)

@inline cross_realm_flux(i, j, grid, time, flux_tuple::Tuple{<:Any, <:Any, <:Any, <:Any}, args...) =
    cross_realm_flux(i, j, grid, time, flux_tuple[1], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[2], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[3], args...) +
    cross_realm_flux(i, j, grid, time, flux_tuple[4], args...)

#####
##### Air-sea differences
#####

@inline air_sea_difference(i, j, grid, time, air, sea) = stateindex(air, i, j, 1, time) -
                                                         stateindex(sea, i, j, 1, time)

struct RelativeUVelocity end
struct RelativeVVelocity end

@inline function air_sea_difference(i, j, grid, time, ::RelativeUVelocity, atmos_state, ocean_state)
    uₐ = atmos_state.u
    uₒ = ocean_state.u
    return air_sea_difference(i, j, grid, time, uₐ, uₒ)
end

@inline function air_sea_difference(i, j, grid, time, ::RelativeVVelocity, atmos_state, ocean_state)
    vₐ = atmos_state.v
    vₒ = ocean_state.v
    return air_sea_difference(i, j, grid, time, vₐ, vₒ)
end

struct SensibleHeat{FT}
    ocean_reference_temperature :: FT
end

SensibleHeat(FT::DataType=Float64; ocean_reference_temperature=273.15) =
    SensibleHeat(convert(FT, ocean_reference_temperature))

@inline function air_sea_difference(i, j, grid, time, Qs::SensibleHeat, atmos_state, ocean_state)
    cₚ = stateindex(atmos_state.cₚ, i, j, 1, time)
    Tₐ = atmos_state.T

    # Compute ocean temperature in degrees K
    Tᵣ = Qs.ocean_reference_temperature 
    Tₒᵢ = stateindex(ocean_state.T, i, j, 1, time)
    Tₒ = Tₒᵢ + Tᵣ
    
    ΔT = air_sea_difference(i, j, grid, time, Tₐ, Tₒ)

    return @inbounds cₚ[i, j, 1] * ΔT
end

struct SpecificHumidity{S}
    saturation_specific_humidity :: S

    @doc """
        SpecificHumidity(FT = Float64;
                               saturation_specific_humidity = LargeYeagerSaturationVaporFraction(FT))

    """
    function SpecificHumidity(FT = Float64;
                              saturation_specific_humidity = LargeYeagerSaturationVaporFraction(FT))
        S = typeof(saturation_specific_humidity)
        return new{S}(saturation_specific_humidity)
    end
end

struct LargeYeagerSaturationVaporFraction{FT}
    q₀ :: FT
    c₁ :: FT
    c₂ :: FT
    reference_temperature :: FT
end

"""
    LargeYeagerSaturationVaporFraction(FT = Float64;
                                       q₀ = 0.98,
                                       c₁ = 640380,
                                       c₂ = -5107.4,
                                       reference_temperature = 273.15)

"""
function LargeYeagerSaturationVaporFraction(FT = Float64;
                                            q₀ = 0.98,
                                            c₁ = 640380,
                                            c₂ = -5107.4,
                                            reference_temperature = 273.15)

    return LargeYeagerSaturationVaporFraction(convert(FT, q₀),
                                              convert(FT, c₁),
                                              convert(FT, c₂),
                                              convert(FT, reference_temperature))
end

@inline function saturation_specific_humidity(i, j, grid, time,
                                           ratio::LargeYeagerSaturationVaporFraction,
                                           atmos_state, ocean_state)

    Tₒ = stateindex(ocean_state.T, i, j, 1, time)
    ρₐ = stateindex(atmos_state.ρ, i, j, 1, time)
    Tᵣ = ratio.reference_temperature
    q₀ = ratio.q₀
    c₁ = ratio.c₁
    c₂ = ratio.c₂

    return q₀ * c₁ * exp(-c₂ / (Tₒ + Tᵣ))
end

@inline function air_sea_difference(i, j, grid, time, diff::SpecificHumidity, atmos_state, ocean_state)
    vapor_fraction = diff.saturation_specific_humidity 
    qₐ = stateindex(atmos_state.q, i, j, 1, time)
    qₛ = saturation_specific_humidity(i, j, grid, time, vapor_fraction, atmos_state, ocean_state)
    return qₐ - qₛ
end

struct LatentHeat{Q, FT}
    specific_humidity_difference :: Q
    vaporization_enthalpy :: FT
end

"""
    LatentHeat(FT = Float64;
               vaporization_enthalpy = 2.5e3 # J / g
               specific_humidity_difference = SpecificHumidity(FT))

"""
function LatentHeat(FT = Float64;
                    vaporization_enthalpy = 2.5e3, # J / g
                    specific_humidity_difference = SpecificHumidity(FT))

    vaporization_enthalpy = convert(FT, vaporization_enthalpy)
    return LatentHeat(specific_humidity_difference, vaporization_enthalpy)
end

@inline function air_sea_difference(i, j, grid, time, diff::LatentHeat, atmos, ocean)
    Δq = air_sea_difference(i, j, grid, time, diff.specific_humidity_difference, atmos, ocean)
    Λᵥ = diff.vaporization_enthalpy
    return Λᵥ * Δq
end

#####
##### Bulk velocity scales
#####

struct RelativeVelocityScale end
# struct AtmosphereOnlyVelocityScale end

@inline function bulk_velocity_scaleᶠᶜᶜ(i, j, grid, time, ::RelativeVelocityScale, atmos_state, ocean_state)
    uₐ = atmos_state.u
    vₐ = atmos_state.v
    uₒ = ocean_state.u
    vₒ = ocean_state.v
    Δu = stateindex(uₐ, i, j, 1, time) - stateindex(uₒ, i, j, 1, time)
    Δv² = ℑxyᶠᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)
    return sqrt(Δu^2 + Δv²)
end

@inline function bulk_velocity_scaleᶜᶠᶜ(i, j, grid, time, ::RelativeVelocityScale, atmos_state, ocean_state)
    uₐ = atmos_state.u
    vₐ = atmos_state.v
    uₒ = ocean_state.u
    vₒ = ocean_state.v
    Δu² = ℑxyᶜᶠᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv = stateindex(vₐ, i, j, 1, time) - stateindex(vₒ, i, j, 1, time)
    return sqrt(Δu² + Δv^2)
end

@inline function bulk_velocity_scaleᶜᶜᶜ(i, j, grid, time, ::RelativeVelocityScale, atmos_state, ocean_state)
    uₐ = atmos_state.u
    vₐ = atmos_state.v
    uₒ = ocean_state.u
    vₒ = ocean_state.v
    Δu² = ℑxᶜᵃᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv² = ℑyᵃᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)
    return sqrt(Δu² + Δv²)
end

#####
##### Convenience containers for surface fluxes
##### 
##### "Cross realm fluxes" can refer to the flux _data_ (ie, fields representing
##### the total flux for a given variable), or to the flux _components_ / formula.
#####

struct CrossRealmFluxes{M, H, T}
    momentum :: M
    heat :: H
    tracers :: T
end

CrossRealmFluxes(; momentum=nothing, heat=nothing, tracers=nothing) =
    CrossRealmFluxes(momentum, heat, tracers)

Base.summary(osf::CrossRealmFluxes) = "CrossRealmFluxes"
Base.show(io::IO, osf::CrossRealmFluxes) = print(io, summary(osf))

=#

