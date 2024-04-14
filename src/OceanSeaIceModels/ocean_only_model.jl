const OceanOnlyModel = OceanSeaIceModel{Nothing}

#####
##### No ice-ocean fluxes in this model!!
#####

function time_step!(coupled_model::OceanOnlyModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean

    # Be paranoid and update state at iteration 0
    coupled_model.clock.iteration == 0 && update_state!(coupled_model, callbacks)

    time_step!(ocean)

    tick!(coupled_model.clock, ocean.Δt) # An Ocean-only model advances with the ocean time-step!
    update_state!(coupled_model, callbacks; compute_tendencies)
    
    return nothing
end

