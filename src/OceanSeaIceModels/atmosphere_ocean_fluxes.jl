function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere

    momentum_fluxes = (u = surface_flux(ocean.model.velocities.u),
                       v = surface_flux(ocean.model.velocities.v))

    tracers = ocean.model.tracers
    tracer_fluxes = NamedTuple(name => surface_flux(tracers[name]) for name in keys(tracers))

    #=
    ε  = coupled_model.ocean_emissivity
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    I₀ = coupled_model.solar_insolation
    =#

    grid = ocean.model.grid
    arch = architecture(grid)
    clock = ocean.model.clock
    ocean_velocities = surface_velocities(ocean)
    atmosphere_velocities = surface_velocities(atmosphere)
    ice_thickness = coupled_model.sea_ice.model.ice_thickness

    #=
    launch!(arch, grid, :xy, _compute_atmosphere_ocean_fluxes!,
            momentum_fluxes, tracer_fluxes,
            ocean_velocities, atmosphere_velocities, ice_thickness)
    =#

    return nothing
end

@kernel function _compute_atmosphere_ocean_fluxes!(momentum_fluxes,
                                                   tracer_fluxes,
                                                   ocean_velocities,
                                                   atmosphere_velocities,
                                                   ice_thickness)

    i, j = @index(Global, NTuple)

    τˣ = momentum_fluxes.u
    τʸ = momentum_fluxes.v

    Cd = 1e-3
    ρₐ = 1.2
    ρₒ = 1020

    @inbounds begin
        ua = atmosphere_velocities.u[i, j, 1]
        va = atmosphere_velocities.v[i, j, 1]

        τˣ[i, j, 1] = ρₐ / ρₒ * Cd * ua * sqrt(ua^2 + va^2)
        τʸ[i, j, 1] = ρₐ / ρₒ * Cd * va * sqrt(ua^2 + va^2)
    end
end

