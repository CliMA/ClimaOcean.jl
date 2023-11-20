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

function extract_top_surface_fluxes(model::HydrostaticFreeSurfaceModel)
    u_flux = surface_flux(model.velocities.u)
    v_flux = surface_flux(model.velocities.v)

    ocean_momentum_fluxes = (u = u_flux, v = v_flux)

    ocean_tracers = model.tracers
    ocean_tracer_fluxes = NamedTuple(name => surface_flux(ocean_tracers[name])
                                     for name in keys(ocean_tracers))

    ocean_fluxes = (momentum = ocean_momentum_fluxes,
                    tracers = ocean_tracer_fluxes)

    return ocean_fluxes
end

#####
##### Total flux across each surface
#####

struct OceanSeaIceSurfaceFluxes{O, IT, IB}
    ocean :: O
    sea_ice_top :: IT
    sea_ice_bottom :: IB
end

function OceanSeaIceSurfaceFluxes(ocean, sea_ice=nothing)
    ocean_fluxes = extract_top_surface_fluxes(ocean.model)
            
    if isnothing(sea_ice)
        sea_ice_top_fluxes = nothing
        sea_ice_bottom_fluxes = nothing
    else
        sea_ice_top_fluxes = extract_top_surface_fluxes(ice.model)
        sea_ice_bottom_fluxes = extract_bottom_surface_fluxes(ice.model)
    end

    return OceanSeaIceSurfaceFluxes(ocean_fluxes,
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

    surfaces = OceanSeaIceSurfaceFluxes(ocean_simulation, sea_ice_simulation)

    return CrossRealmFluxes(surfaces,
                            radiation,
                            atmosphere_ocean,
                            atmosphere_sea_ice,
                            sea_ice_ocean)
end


