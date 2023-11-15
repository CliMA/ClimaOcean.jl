function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere

    momentum_fluxes = (u = surface_flux(ocean.model.velocities.u),
                       v = surface_flux(ocean.model.velocities.v))

    tracers = ocean.model.tracers
    tracer_fluxes = NamedTuple(name => surface_flux(tracers[name]) for name in keys(tracers))

    grid = ocean.model.grid
    arch = architecture(grid)
    clock = ocean.model.clock
    ocean_velocities = surface_velocities(ocean)
    ocean_reference_density = coupled_model.ocean_reference_density
    atmosphere_velocities = surface_velocities(atmosphere)
    ice_thickness = coupled_model.sea_ice.model.ice_thickness

    launch!(arch, grid, :xy, _compute_atmosphere_ocean_fluxes!,
            grid, clock, momentum_fluxes, tracer_fluxes,
            ocean_velocities, atmosphere_velocities, ocean_reference_density,
            ice_thickness)

    return nothing
end


@kernel function _compute_atmosphere_ocean_fluxes!(grid,
                                                   clock,
                                                   momentum_fluxes,
                                                   tracer_fluxes,
                                                   ocean_velocities,
                                                   atmosphere_velocities,
                                                   ocean_reference_density,
                                                   ice_thickness)

    i, j = @index(Global, NTuple)

    τˣ = momentum_fluxes.u
    τʸ = momentum_fluxes.v
    time = Time(clock.time)

    Uₒ = ocean_velocities
    Uₐ = atmosphere_velocities
    cᴰ = 1e-3
    ρₐ = 1.2
    ρₒ = ocean_reference_density

    @inbounds begin
        τˣ[i, j, 1] = x_atmosphere_ocean_momentum_flux(i, j, grid, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
        τʸ[i, j, 1] = y_atmosphere_ocean_momentum_flux(i, j, grid, cᴰ, ρₒ, ρₐ, Uₒ, Uₐ)
    end
end

