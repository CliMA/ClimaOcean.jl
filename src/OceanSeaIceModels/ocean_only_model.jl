const OceanOnlyModel = OceanSeaIceModel{Nothing}
const OceanCappedSeaIceModel = OceanSeaIceModel{MinimumTemperatureSeaIce}

#####
##### No ice-ocean fluxes in this models!!
#####

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: compute_sea_ice_ocean_fluxes!

compute_sea_ice_ocean_fluxes!(::OceanOnlyModel) = nothing
compute_sea_ice_ocean_fluxes!(::OceanCappedSeaIceModel) = nothing

function time_step!(coupled_model::Union{OceanOnlyModel, OceanCappedSeaIceModel}, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean

    # Be paranoid and update state at iteration 0
    coupled_model.clock.iteration == 0 && update_state!(coupled_model, callbacks)

    time_step!(ocean)

    tick!(coupled_model.clock, ocean.Δt) # An Ocean-only model advances with the ocean time-step!
    update_state!(coupled_model, callbacks; compute_tendencies)
    
    return nothing
end

function update_state!(coupled_model::Union{OceanOnlyModel, OceanCappedSeaIceModel}, callbacks=[]; compute_tendencies=false)
    time = Time(coupled_model.clock.time)
    update_model_field_time_series!(coupled_model.atmosphere, time)
    compute_atmosphere_ocean_fluxes!(coupled_model) 
    return nothing
end

