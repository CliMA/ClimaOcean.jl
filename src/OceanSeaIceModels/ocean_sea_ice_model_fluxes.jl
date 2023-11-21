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
##### Extractors for differnet models (maybe this belongs in the model repo's)
#####

function extract_top_surface_fluxes(model::HydrostaticFreeSurfaceModel)
    u_flux = surface_flux(model.velocities.u)
    v_flux = surface_flux(model.velocities.v)

    ocean_momentum_fluxes = (u = u_flux, v = v_flux)

    ocean_tracers = model.tracers

    ocean_tracer_fluxes = NamedTuple(name => surface_flux(ocean_tracers[name])
                                     for name in keys(ocean_tracers)
                                     if surface_flux(ocean_tracers[name]) isa AbstractArray)

    ocean_fluxes = CrossRealmFluxes(momentum = ocean_momentum_fluxes,
                                    tracers = ocean_tracer_fluxes)

    return ocean_fluxes
end

extract_top_surface_fluxes(model::SlabSeaIceModel) = nothing
extract_bottom_surface_fluxes(model::SlabSeaIceModel) = nothing

#####
##### Total flux across each surface
#####

struct OceanSeaIceSurfaces{O, IT, IB}
    ocean :: O
    sea_ice_top :: IT
    sea_ice_bottom :: IB
end

Base.summary(osis::OceanSeaIceSurfaces) = "OceanSeaIceSurfaces"
Base.show(io::IO, osis::OceanSeaIceSurfaces) = print(io, summary(osis))

function OceanSeaIceSurfaces(ocean, sea_ice=nothing)
    ocean_fluxes = extract_top_surface_fluxes(ocean.model)
            
    if isnothing(sea_ice)
        sea_ice_top_fluxes = nothing
        sea_ice_bottom_fluxes = nothing
    else
        sea_ice_top_fluxes = extract_top_surface_fluxes(sea_ice.model)
        sea_ice_bottom_fluxes = extract_bottom_surface_fluxes(sea_ice.model)
    end

    return OceanSeaIceSurfaces(ocean_fluxes,
                               sea_ice_top_fluxes,  
                               sea_ice_bottom_fluxes)
end

#####
##### Container for organizing information related to fluxes
#####

struct OceanSeaIceModelFluxes{S, R, AO, AI, IO}
    surfaces :: S
    radiation :: R
    atmosphere_ocean :: AO
    atmosphere_sea_ice :: AI
    sea_ice_ocean :: IO
end

function OceanSeaIceModelFluxes(ocean, sea_ice=nothing, atmosphere=nothing;
                                radiation = nothing,
                                atmosphere_ocean = nothing,
                                atmosphere_sea_ice = nothing,
                                sea_ice_ocean = nothing)

    surfaces = OceanSeaIceSurfaces(ocean, sea_ice)

    if isnothing(atmosphere_ocean) # defaults
        FT = eltype(ocean.model.grid)
        τˣ = BulkFormula(FT, transfer_coefficient=1e-3)
        τʸ = BulkFormula(FT, transfer_coefficient=1e-3)
        momentum_flux_formulae = (u=τˣ, v=τʸ)
        atmosphere_ocean = CrossRealmFluxes(momentum = momentum_flux_formulae)
    end

    return OceanSeaIceModelFluxes(surfaces,
                                  radiation,
                                  atmosphere_ocean,
                                  atmosphere_sea_ice,
                                  sea_ice_ocean)
end

Base.summary(crf::OceanSeaIceModelFluxes) = "OceanSeaIceModelFluxes"
Base.show(io::IO, crf::OceanSeaIceModelFluxes) = print(io, summary(crf))

#####
##### Bulk formula
#####

struct RelativeVelocityScale end
struct AtmosphereOnlyVelocityVelocityScale end

struct BulkFormula{T, CD}
    velocity_scale :: T
    transfer_coefficient :: CD
end

function BulkFormula(FT=Float64;
                     velocity_scale = RelativeVelocityScale(),
                     transfer_coefficient = 1e-3)

    return BulkFormula(velocity_scale,
                       convert(FT, transfer_coefficient))
end

