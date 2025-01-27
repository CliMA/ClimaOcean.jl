using .InterfaceComputations:
    compute_atmosphere_ocean_fluxes!,
    compute_sea_ice_ocean_fluxes!,
    compute_net_ocean_fluxes!,
    interpolate_atmospheric_state!

using ClimaSeaIce: SeaIceModel, SeaIceThermodynamics

using Printf

function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice

    # Be paranoid and update state at iteration 0
    coupled_model.clock.iteration == 0 && update_state!(coupled_model, callbacks)

    # Eventually, split out into OceanOnlyModel
    if sea_ice isa SeaIceSimulation
        h = sea_ice.model.ice_thickness
        fill_halo_regions!(h)

        # Initialization
        if coupled_model.clock.iteration == 0
            @info "Initializing coupled model ice thickness..."
            h⁻ = coupled_model.interfaces.sea_ice_ocean_interface.previous_ice_thickness
            hⁿ = coupled_model.sea_ice.model.ice_thickness
            parent(h⁻) .= parent(hⁿ)
        end

        sea_ice.Δt = Δt
        thermodynamic_sea_ice_time_step!(coupled_model)
    end

    # TODO after ice time-step:
    #  - Adjust ocean heat flux if the ice completely melts?
    ocean.Δt = Δt
    time_step!(ocean)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    tick!(coupled_model.clock, Δt)
    update_state!(coupled_model, callbacks; compute_tendencies)

    return nothing
end

function update_state!(coupled_model::OceanSeaIceModel, callbacks=[]; compute_tendencies=true)
    time = Time(coupled_model.clock.time)
    update_model_field_time_series!(coupled_model.atmosphere, time)
    interpolate_atmospheric_state!(coupled_model)

    # Compute interface states
    compute_atmosphere_ocean_fluxes!(coupled_model)
    compute_atmosphere_sea_ice_fluxes!(coupled_model)
    compute_sea_ice_ocean_fluxes!(coupled_model)

    # compute_net_atmosphere_fluxes!(coupled_model)
    compute_net_ocean_fluxes!(coupled_model)
    #compute_net_sea_ice_fluxes!(coupled_model)

    return nothing
end

function thermodynamic_sea_ice_time_step!(coupled_model)
    sea_ice = coupled_model.sea_ice
    model = sea_ice.model
    Δt = sea_ice.Δt
    grid = coupled_model.ocean.model.grid
    arch = architecture(grid)
    clock = model.clock
    thermodynamics = model.ice_thermodynamics
    ice_thickness = model.ice_thickness
    ice_concentration = model.ice_concentration
    top_external_heat_flux = model.external_heat_fluxes.top
    bottom_external_heat_flux = model.external_heat_fluxes.bottom
    ocean_salinity = coupled_model.ocean.model.tracers.S

    launch!(arch, grid, :xy, update_thickness!,
            ice_thickness,
            grid, Δt,
            ice_concentration,
            ocean_salinity,
            thermodynamics,
            top_external_heat_flux,
            bottom_external_heat_flux,
            clock)

    return nothing
end

@kernel function update_thickness!(ice_thickness,
                                   grid, Δt,
                                   ice_concentration,
                                   ocean_salinity,
                                   thermodynamics,
                                   top_external_heat_flux,
                                   bottom_external_heat_flux,
                                   clock)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)

    phase_transitions = thermodynamics.phase_transitions
    top_heat_bc = thermodynamics.heat_boundary_conditions.top
    bottom_heat_bc = thermodynamics.heat_boundary_conditions.bottom
    liquidus = phase_transitions.liquidus

    Qi = thermodynamics.internal_heat_flux
    Qu = top_external_heat_flux
    Qb = bottom_external_heat_flux
    Tu = thermodynamics.top_surface_temperature

    @inbounds begin
        hᶜ = thermodynamics.ice_consolidation_thickness
        hᵢ = ice_thickness[i, j, 1]
        ℵᵢ = ice_concentration[i, j, 1]
    end

    # Consolidation criteria
    @inbounds Tuᵢ = Tu[i, j, 1]

    # Bottom temperature at the melting temperature
    Sₒ = @inbounds ocean_salinity[i, j, kᴺ]
    Tbᵢ = SeaIceThermodynamics.melting_temperature(liquidus, Sₒ)
    ℰb = SeaIceThermodynamics.latent_heat(phase_transitions, Tbᵢ)
    ℰu = SeaIceThermodynamics.latent_heat(phase_transitions, Tuᵢ)

    # Retrieve fluxes
    @inbounds begin
        Quᵢ = Qu[i, j, 1]
        Qbᵢ = Qb[i, j, 1]
    end

    # If ice is consolidated, compute tendency for an ice slab; otherwise
    # just add ocean fluxes from frazil ice formation or melting
    # wb = - Qbᵢ / ℰb

    # Clip thickness for thermodynamic computations
    #hᵢ = max(hᶜ, hᵢ)
    k = Qi.parameters.flux.conductivity
    Qiᵢ = - k * (Tuᵢ - Tbᵢ) / hᵢ * (hᵢ > hᶜ) # getflux(Qi, i, j, grid, Tuᵢ, clock, model_fields)

    # Upper (top) and bottom interface velocities
    w_melting  = (Quᵢ - Qiᵢ) / ℰu # < 0 => melting
    w_freezing = +Qiᵢ / ℰb # < 0 => freezing
    w_frazil   = -Qbᵢ / ℰb # < 0 => freezing

    Δh_melting  = w_melting  * Δt * ℵᵢ
    Δh_freezing = w_freezing * Δt * ℵᵢ
    Δh_frazil   = w_frazil   * Δt # frazil flux contributes from entire cell

    if ℵᵢ < 1 && Δh_frazil > 0 # Add ice volume laterally
        # ΔV = h * Δℵ
        ℵ⁺ = ℵᵢ + Δh_frazil / hᶜ
        ℵ⁺ = min(one(ℵ⁺), ℵ⁺)
        Δh_frazil -= ℵ⁺ * hᵢ
    else
        ℵ⁺ = ℵᵢ
    end

    Δh = Δh_frazil + Δh_freezing + Δh_melting

    if Δh > 0
        h⁺ = hᵢ + Δh / ℵ⁺
        h⁺ = max(zero(h⁺), h⁺)
    else
        h⁺ = hᵢ
    end

    # TODO: incorporate minimum_thickness

    @inbounds begin
        ice_thickness[i, j, 1] = h⁺
        ice_concentration[i, j, 1] = ℵ⁺
    end

    #=
    @printf("Tu: %.1f, Δh freeze: %.1e, Δh melt : %.1e, Δh frazil: %.1e \n",
            Tuᵢ, Δh_freezing, Δh_melting, Δh_frazil)

    @printf("h: %.1e, ℵ: %.1e Qu: %.1e, Qi: %.1e, Qb: %.1e \n",
            h⁺, ℵ⁺, Quᵢ, Qiᵢ, Qbᵢ)
    =#

end

