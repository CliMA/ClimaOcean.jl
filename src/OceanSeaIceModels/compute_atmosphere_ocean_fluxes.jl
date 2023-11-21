function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere

    tracers = ocean.model.tracers
    tracer_fluxes = NamedTuple(name => surface_flux(tracers[name]) for name in keys(tracers))

    grid = ocean.model.grid
    arch = architecture(grid)
    clock = ocean.model.clock
    ocean_velocities = surface_velocities(ocean)
    ocean_reference_density = coupled_model.ocean_reference_density
    atmosphere_velocities = surface_velocities(atmosphere)
    ice_thickness = coupled_model.sea_ice.model.ice_thickness

    momentum_flux_fields = coupled_model.fluxes.surfaces.ocean.momentum
    tracer_flux_fields = coupled_model.fluxes.surfaces.ocean.tracers
    momentum_flux_contributions = coupled_model.fluxes.atmosphere_ocean.momentum
    heat_flux_contributions = coupled_model.fluxes.atmosphere_ocean.heat
    tracer_flux_contributions = coupled_model.fluxes.atmosphere_ocean.heat

    atmosphere_ocean_parameters = (
        ρₐ = 1.2,
        ρₒ = 1024,,
        cₚ = 3991.0,
    )

    launch!(arch, grid, :xy, _compute_atmosphere_ocean_fluxes!,
            grid, clock,
            momentum_flux_fields,
            tracer_flux_fields,
            momentum_flux_contributions,
            heat_flux_contributions
            ocean_velocities,
            atmosphere_velocities,
            atmosphere_ocean_parameters,
            ice_thickness)

    return nothing
end

@inline function air_sea_difference(i, j, grid, time, ::RelativeAtmosphereOceanVelocity,
                                    air::SomeKindOfFieldTimeSeries, sea::AbstractArray)

    δ = @inbounds air[i, j, 1, time] - sea[i, j, 1]
    #@show air[i, j, 1, time] time δ
    return δ
end

@kernel function _compute_atmosphere_ocean_fluxes!(grid,
                                                   clock,
                                                   momentum_flux_fields,
                                                   tracer_flux_fields,
                                                   momentum_flux_contributions,
                                                   heat_flux_contributions,
                                                   tracer_flux_contributions,
                                                   ocean_velocities,
                                                   atmosphere_velocities,
                                                   radiation,
                                                   atmosphere_ocean_parameters,
                                                   ice_thickness)

    i, j = @index(Global, NTuple)

    τˣ = momentum_flux_fields.u
    τʸ = momentum_flux_fields.v
    # Q = tracer_fluxes.T
    # F = tracer_fluxes.S

    time = Time(clock.time)

    Uₒ = ocean_velocities
    Uₐ = atmosphere_velocities
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    ρₐ = atmosphere_ocean_parameters.ρₐ
    ρₒ = atmosphere_ocean_parameters.ρₒ
    cₚ = atmosphere_ocean_parameters.cₚ
    time = Time(clock.time)

    u_formula = momentum_flux_contributions.u.transfer_velocity
    v_formula = momentum_flux_contributions.v.transfer_velocity

    cᴰ = momentum_flux_contributions.u.transfer_coefficient

    # Compute transfer velocity scale
    Vᶠᶜᶜ = transfer_velocityᶠᶜᶜ(i, j, grid, time, u_formula, Uₐ, Uₒ)
    Vᶜᶠᶜ = transfer_velocityᶜᶠᶜ(i, j, grid, time, v_formula, Uₐ, Uₒ)

    Δu = air_sea_difference(i, j, grid, time, u_formula, uₐ, uₒ)
    Δv = air_sea_difference(i, j, grid, time, v_formula, vₐ, vₒ)

    @inbounds begin
        τˣ[i, j, 1] = - ρₐ / ρₒ * cᴰ * Δu * Vᶠᶜᶜ
        τʸ[i, j, 1] = - ρₐ / ρₒ * cᴰ * Δv * Vᶜᶠᶜ

        # @show τˣ[i, j, 1] 
        # @show τʸ[i, j, 1]

        Q[i, j, 1] = radiative_fluxes(i, j, grid, time, radiation)
    end
end

@inline function radiative_fluxes(i, j, grid, time, radiation)
    Qˢʷ = radiation.downwelling_shortwave_radiation
    Qˡʷ = radiation.downwelling_longwave_radiation
    α = radiation.ocean_albedo
    return @inbounds (1 - α) * Qˢʷ[i, j, 1, time] + Qˡʷ[i, j, 1, time]
end

