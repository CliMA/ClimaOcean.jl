using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!

const OceanOnlyModel = OceanSeaIceModel{Nothing}
const OceanSimplifiedSeaIceModel = OceanSeaIceModel{<:MinimumTemperatureSeaIce}

const NoSeaIceModel = Union{OceanOnlyModel, OceanSimplifiedSeaIceModel}

#####
##### No ice-ocean fluxes in this models!!
#####

import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: compute_sea_ice_ocean_fluxes!

compute_sea_ice_ocean_fluxes!(::NoSeaIceModel) = nothing

function time_step!(coupled_model::NoSeaIceModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean

    # Be paranoid and update state at iteration 0
    coupled_model.clock.iteration == 0 && update_state!(coupled_model, callbacks)

    time_step!(ocean)

    tick!(coupled_model.clock, ocean.Δt) # An Ocean-only model advances with the ocean time-step!
    update_state!(coupled_model, callbacks; compute_tendencies)
    
    return nothing
end

function update_state!(coupled_model::NoSeaIceModel, callbacks=[]; compute_tendencies=false)
    time = Time(coupled_model.clock.time)
    update_model_field_time_series!(coupled_model.atmosphere, time)
    
    ocean_model = coupled_model.ocean.model

    # Do we really have to do this?
    if !isempty(ocean_model.forcing)
        field_time_series = extract_field_time_series(ocean_model.forcing)
        update_field_time_series!(field_time_series, time)
    end

    compute_atmosphere_ocean_fluxes!(coupled_model) 
    return nothing
end

