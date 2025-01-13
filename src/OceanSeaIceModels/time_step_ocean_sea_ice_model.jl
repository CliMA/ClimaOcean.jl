using .CrossRealmFluxes: compute_atmosphere_ocean_fluxes!, compute_sea_ice_ocean_fluxes!, interpolate_atmospheric_state!

using ClimaSeaIce: SeaIceModel

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
            h⁻ = coupled_model.fluxes.previous_ice_thickness
            hⁿ = coupled_model.sea_ice.model.ice_thickness
            parent(h⁻) .= parent(hⁿ)
        end

        sea_ice.Δt = Δt
        time_step!(sea_ice)
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
    compute_atmosphere_ocean_fluxes!(coupled_model)

    # The atmospheric state is interpolated on the ocean grid in the previous 
    # step (compute_atmosphere_ocean_fluxes!). We reuse the same interpolated
    # state for the sea ice model assuming that the ocean and the sea ice live on 
    # the same horizontal grid.
    compute_atmosphere_sea_ice_fluxes!(coupled_model)
    compute_sea_ice_ocean_fluxes!(coupled_model)

    return nothing
end
