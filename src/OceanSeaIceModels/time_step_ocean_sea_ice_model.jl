using .InterfaceComputations:
    compute_atmosphere_ocean_fluxes!,
    compute_sea_ice_ocean_fluxes!

using ClimaSeaIce: SeaIceModel, SeaIceThermodynamics
using Oceananigans.Grids: φnode
using Printf

function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    atmosphere = coupled_model.atmosphere

    # Eventually, split out into OceanOnlyModel
    !isnothing(sea_ice) && time_step!(sea_ice, Δt)
    
    # TODO after ice time-step:
    #  - Adjust ocean heat flux if the ice completely melts?
    !isnothing(ocean) && time_step!(ocean, Δt)

    # Time step the atmosphere
    !isnothing(atmosphere) && time_step!(atmosphere, Δt)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    tick!(coupled_model.clock, Δt)
    update_state!(coupled_model, callbacks; compute_tendencies)

    return nothing
end

function update_state!(coupled_model::OceanSeaIceModel, callbacks=[]; compute_tendencies=true)

    # The three components
    ocean      = coupled_model.ocean
    sea_ice    = coupled_model.sea_ice
    atmosphere = coupled_model.atmosphere

    exchanger = coupled_model.interfaces.exchanger
    grid      = exchanger.grid
    
    # This function needs to be specialized to allow different component models
    interpolate_state!(exchanger.atmosphere, grid, atmosphere, coupled_model)
    interpolate_state!(exchanger.ocean,      grid, ocean,      coupled_model)
    interpolate_state!(exchanger.sea_ice,    grid, sea_ice,    coupled_model)

    # Compute interface states
    compute_atmosphere_ocean_fluxes!(coupled_model)
    compute_atmosphere_sea_ice_fluxes!(coupled_model)
    compute_sea_ice_ocean_fluxes!(coupled_model)

    # This function needs to be specialized to allow different component models
    update_net_fluxes!(coupled_model, atmosphere)
    update_net_fluxes!(coupled_model, ocean)
    update_net_fluxes!(coupled_model, sea_ice)

    return nothing
end
