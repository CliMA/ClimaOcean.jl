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


