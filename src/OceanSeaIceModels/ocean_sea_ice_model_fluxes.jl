using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
using ClimaSeaIce.SlabSeaIceModels: SlabSeaIceModel

#####
##### Utilities
#####

function surface_flux(f::Field)
    top_bc = f.boundary_conditions.top
    if top_bc isa BoundaryCondition{<:Oceananigans.BoundaryConditions.Flux}
        return top_bc.condition
    else
        return nothing
    end
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

#####
##### Container for organizing information related to fluxes
#####

struct RelativeVelocityScale end
struct AtmosphereOnlyVelocityScale end

struct OceanSeaIceModelFluxes{U, R, AO, ASI, SIO}
    bulk_velocity_scale :: U
    downwelling_radiation :: R
    atmosphere_ocean :: AO
    atmosphere_sea_ice :: ASI
    sea_ice_ocean :: SIO
end

function OceanSeaIceModelFluxes(FT=Float64;
                                bulk_velocity_scale = RelativeVelocityScale(),
                                downwelling_radiation = nothing,
                                atmosphere_ocean = nothing,
                                atmosphere_sea_ice = nothing,
                                sea_ice_ocean = nothing)

    if isnothing(atmosphere_ocean) # defaults
        τˣ = BulkFormula(RelativeVelocity(), 1e-3)
        τʸ = BulkFormula(RelativeVelocity(), 1e-3)
        momentum_flux_formulae = (u=τˣ, v=τʸ)

        evaporation = BulkFormula(SpecificHumidity, 1e-3)
                

        atmosphere_ocean = CrossRealmFluxes(momentum = momentum_flux_formulae)
    end

    return OceanSeaIceModelFluxes(bulk_velocity_scale,
                                  downwelling_radiation,
                                  atmosphere_ocean,
                                  atmosphere_sea_ice,
                                  sea_ice_ocean)
end

Base.summary(crf::OceanSeaIceModelFluxes) = "OceanSeaIceModelFluxes"
Base.show(io::IO, crf::OceanSeaIceModelFluxes) = print(io, summary(crf))

#####
##### Bulk formula
#####

"""
    BulkFormula(air_sea_difference, transfer_coefficient)

The basic structure of a flux `J` computed by a bulk formula is:

```math
J = ρₐ * C * Δc * ΔU
```

where `ρₐ` is the density of air, `C` is the `transfer_coefficient`,
`Δc` is the air_sea_difference, and `ΔU` is the bulk velocity scale.
"""
struct BulkFormula{F, CD}
    air_sea_difference :: F
    transfer_coefficient :: CD
end

@inline function tracer_flux(i, j, grid, time, formula::BulkFormula, ΔU, atmosphere_state, ocean_state)
    ρₐ = atmosphere_state.density
    C = formula.transfer_coefficient
    Δc = air_sea_difference(i, j, grid, time, formula.air_sea_difference, atmosphere_state, ocean_state)
    return ρₐ * C * Δc * ΔU
end

@inline tracer_flux(i, j, grid, time, flux::NamedTuple, args...) =
    tracer_flux(i, j, grid, time, values(flux))

@inline tracer_flux(i, j, grid, time, flux::Tuple{<:Any, <:Any}, args...) =
    tracer_flux(i, j, grid, time, flux[1], args...) +
    tracer_flux(i, j, grid, time, flux[2], args...)

@inline tracer_flux(i, j, grid, time, fts::SKOFTS, args...) =
    @inbounds fts[i, j, 1, time]

#####
##### Air-sea differences
#####

struct RelativeVelocity end

struct MassSpecificHumidity{S}
    saturation :: S
end

struct LargeYeagerSaturation{FT}
    c1 :: FT
    c2 :: FT
    reference_temperature:: FT
end

function LargeYeagerSaturation(FT=Float64;
                               c1 = 0.98 * 640380,
                               c2 = -5107.4,
                               reference_temperature = 273.15)
    return LargeYeagerSaturation(convert(FT, c1),
                                 convert(FT, c2),
                                 convert(FT, reference_temperature))
end

#=
MassSpecificHumidity(FT=Float64; saturation = LargeYeagerSaturation(FT)) =
    MassSpecificHumidity(saturation)
                              
struct InternalEnergy{C}
    atmosphere_specific_heat :: C
end

InternalEnergy(FT::DataType; atmosphere_specific_heat=1000.5) =
    InternalEnergy(convert(FT,  atmosphere_specific_heat))
=#

@inline function air_sea_difference(i, j, grid, time, ::MassSpecificHumidity, atmos, ocean)
    return air_sea_difference(i, j, grid, time, atmos, ocean)
end

@inline air_sea_difference(i, j, grid, time, air::SKOFTS, sea::AbstractArray) =
    @inbounds air[i, j, 1, time] - sea[i, j, 1]

#####
##### Bulk velocity scales
#####

@inline function bulk_velocity_scaleᶠᶜᶜ(i, j, grid, time, ::RelativeVelocityScale, Uₐ, Uₒ)
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    Δu = @inbounds uₐ[i, j, 1, time] - uₒ[i, j, 1]
    Δv² = ℑxyᶠᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)

    return sqrt(Δu^2 + Δv²)
end

@inline function bulk_velocity_scaleᶜᶠᶜ(i, j, grid, time, ::RelativeVelocityScale, Uₐ, Uₒ)
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    Δu² = ℑxyᶜᶠᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv = @inbounds vₐ[i, j, 1, time] - vₒ[i, j, 1]

    return sqrt(Δu² + Δv^2)
end

@inline function bulk_velocity_scaleᶜᶜᶜ(i, j, grid, time, ::RelativeVelocityScale, Uₐ, Uₒ)
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    Δu² = ℑxᶜᵃᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv² = ℑyᵃᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)

    return sqrt(Δu² + Δv²)
end

