using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
using ClimaSeaIce.SlabSeaIceModels: SlabSeaIceModel

#####
##### Utilities
#####

function surface_velocities(ocean::Simulation{<:HydrostaticFreeSurfaceModel})
    grid = ocean.model.grid
    Nz = size(grid, 3)
    u = interior(ocean.model.velocities.u, :, :, Nz)
    v = interior(ocean.model.velocities.v, :, :, Nz)
    w = interior(ocean.model.velocities.w, :, :, Nz+1)
    return (; u, v, w)
end

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

struct OceanSurfaceFluxes{M, T}
    momentum :: M
    tracers :: T
end

Base.summary(osf::OceanSurfaceFluxes) = "OceanSurfaceFluxes"
Base.show(io::IO, osf::OceanSurfaceFluxes) = print(io, summary(osf))

struct SeaIceSurfaceFluxes{M, T}
    momentum :: M
    tracers :: T
end

Base.summary(sisf::SeaIceSurfaceFluxes) = "SeaIceSurfaceFluxes"
Base.show(io::IO, sisf::SeaIceSurfaceFluxes) = print(io, summary(sisf))

function extract_top_surface_fluxes(model::HydrostaticFreeSurfaceModel)
    u_flux = surface_flux(model.velocities.u)
    v_flux = surface_flux(model.velocities.v)

    ocean_momentum_fluxes = (u = u_flux, v = v_flux)

    ocean_tracers = model.tracers

    ocean_tracer_fluxes = NamedTuple(name => surface_flux(ocean_tracers[name])
                                     for name in keys(ocean_tracers)
                                     if surface_flux(ocean_tracers[name]) isa AbstractArray)

    ocean_fluxes = OceanSurfaceFluxes(ocean_momentum_fluxes,
                                      ocean_tracer_fluxes)

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
##### Cross-realm fluxes
#####

struct CrossRealmFluxes{S, R, AO, AI, IO}
    surfaces :: S
    radiation :: R
    atmosphere_ocean :: AO
    atmosphere_sea_ice :: AI
    sea_ice_ocean :: IO
end

function CrossRealmFluxes(ocean_simulation, sea_ice_simulation=nothing;
                          radiation = nothing,
                          atmosphere_ocean = nothing,
                          atmosphere_sea_ice = nothing,
                          sea_ice_ocean = nothing)

    surfaces = OceanSeaIceSurfaces(ocean_simulation, sea_ice_simulation)

    return CrossRealmFluxes(surfaces,
                            radiation,
                            atmosphere_ocean,
                            atmosphere_sea_ice,
                            sea_ice_ocean)
end

Base.summary(crf::CrossRealmFluxes) = "CrossRealmFluxes"
Base.show(io::IO, crf::CrossRealmFluxes) = print(io, summary(crf))

#####
##### CrossRealmFlux
#####

struct RelativeAtmosphereOceanVelocity end
struct AtmosphereVelocity end

#####
##### Bulk formula
#####

struct BulkFormula{T, CD}
    transfer_velocity :: T
    transfer_coefficient :: CD
end

function BulkFormula(FT=Float64;
                     transfer_velocity = RelativeAtmosphereOceanVelocity(),
                     transfer_coefficient = convert(FT, 1e-3))

    return BulkFormula(transfer_velocity, transfer_coefficient)
end

#####
##### Abstraction for fluxes across the realms
#####

struct CrossRealmFlux{EQ, F}
    formula :: EQ
    flux :: F
end

"""
    CrossRealmFlux(flux_field; formula = nothing)

May the realms communicate.
"""
function CrossRealmFlux(flux_field; formula = nothing)

    if isnothing(formula) # constant coefficient then
        formula = BulkFormula(eltype(flux_field))
    end
                        
    return CrossRealmFlux(formula, flux_field)
end

