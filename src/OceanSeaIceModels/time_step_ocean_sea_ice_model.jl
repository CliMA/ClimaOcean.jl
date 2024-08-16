using .CrossRealmFluxes: compute_atmosphere_ocean_fluxes!, compute_sea_ice_ocean_fluxes!
function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice

    # Be paranoid and update state at iteration 0
    coupled_model.clock.iteration == 0 && update_state!(coupled_model, callbacks)

    ocean.Δt = Δt

    # TODO after ice time-step:
    #   - Adjust ocean heat flux if the ice completely melts?

    time_step!(ocean)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    tick!(coupled_model.clock, Δt)
    
    update_state!(coupled_model, callbacks; compute_tendencies)
    
    return nothing
end

function update_state!(coupled_model::OceanSeaIceModel, callbacks=[]; compute_tendencies=false)
    time = Time(coupled_model.clock.time)
    update_model_field_time_series!(coupled_model.atmosphere, time)
    compute_atmosphere_ocean_fluxes!(coupled_model) 
    compute_sea_ice_ocean_fluxes!(coupled_model)
    #compute_atmosphere_sea_ice_fluxes!(coupled_model)
    return nothing
end

