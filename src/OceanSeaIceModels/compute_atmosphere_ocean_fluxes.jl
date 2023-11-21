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
    momentum_flux_formulae = coupled_model.fluxes.atmosphere_ocean.momentum

    launch!(arch, grid, :xy, _compute_atmosphere_ocean_fluxes!,
            grid, clock,
            momentum_flux_fields,
            momentum_flux_formulae,
            ocean_velocities,
            atmosphere_velocities,
            ice_thickness)

    return nothing
end

@inline function air_sea_difference(i, j, grid, time, ::RelativeAtmosphereOceanVelocity,
                                    air::SomeKindOfFieldTimeSeries, sea::AbstractArray)

    δ = @inbounds air[i, j, 1, time] - sea[i, j, 1]
    @show air[i, j, 1, time] time δ
    return δ
end

@kernel function _compute_atmosphere_ocean_fluxes!(grid,
                                                   clock,
                                                   momentum_flux_fields,
                                                   momentum_flux_formula,
                                                   ocean_velocities,
                                                   atmosphere_velocities,
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

    cᴰ = 1e-3
    ρₐ = 1.2
    ρₒ = 1024
    time = Time(clock.time)

    u_formula = momentum_flux_formula.u.transfer_velocity
    v_formula = momentum_flux_formula.v.transfer_velocity

    cᴰ = momentum_flux_formula.u.transfer_coefficient

    # Compute transfer velocity scale
    Vᶠᶜᶜ = transfer_velocityᶠᶜᶜ(i, j, grid, time, u_formula, Uₐ, Uₒ)
    Vᶜᶠᶜ = transfer_velocityᶜᶠᶜ(i, j, grid, time, v_formula, Uₐ, Uₒ)

    Δu = air_sea_difference(i, j, grid, time, u_formula, uₐ, uₒ)
    Δv = air_sea_difference(i, j, grid, time, v_formula, vₐ, vₒ)

    @show uₒ[i, j, 1]
    @show vₒ[i, j, 1]
    @show Vᶠᶜᶜ
    @show Vᶜᶠᶜ
    @show Δu
    @show Δv

    @inbounds begin
        τˣ[i, j, 1] = - ρₐ / ρₒ * cᴰ * Δu * Vᶠᶜᶜ
        τʸ[i, j, 1] = - ρₐ / ρₒ * cᴰ * Δv * Vᶜᶠᶜ

        @show τˣ[i, j, 1] 
        @show τʸ[i, j, 1]

        # Q[i, j, 1] = ρₐ
    end
end

