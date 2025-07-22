using .InterfaceComputations:
    compute_atmosphere_ocean_fluxes!,
    compute_sea_ice_ocean_fluxes!,
    compute_net_ocean_fluxes!,
    compute_net_sea_ice_fluxes!,
    interpolate_atmosphere_state!

using ClimaSeaIce: SeaIceModel, SeaIceThermodynamics
using Oceananigans.Grids: φnode
using Oceananigans.Simulations: TimeStepWizard

using Printf

function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    atmosphere = coupled_model.atmosphere

    # Eventually, split out into OceanOnlyModel
    !isnothing(sea_ice) && time_step!(sea_ice, Δt)
    
    # TODO after ice time-step:
    #  - Adjust ocean heat flux if the ice completely melts?
    time_step!(ocean, Δt)

    # Time step the atmosphere
    time_step!(atmosphere, Δt)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    tick!(coupled_model.clock, Δt)
    update_state!(coupled_model, callbacks; compute_tendencies)

    return nothing
end

function update_state!(coupled_model::OceanSeaIceModel, callbacks=[]; compute_tendencies=true)

    # This function needs to be specialized to allow different atmospheric models
    interpolate_atmosphere_state!(coupled_model.interfaces, coupled_model.atmosphere, coupled_model)

    # Compute interface states
    compute_atmosphere_ocean_fluxes!(coupled_model)
    compute_atmosphere_sea_ice_fluxes!(coupled_model)
    compute_sea_ice_ocean_fluxes!(coupled_model)

    # This function needs to be specialized to allow different atmospheric models
    compute_net_atmosphere_fluxes!(coupled_model)
    compute_net_ocean_fluxes!(coupled_model)
    compute_net_sea_ice_fluxes!(coupled_model)

    return nothing
end

function (wizard::TimeStepWizard)(simulation::Simulation{<:OceanSeaIceModel}) 
    model = simulation.model
    ocean_Δt = wizard(model.ocean)
    
    sea_ice_Δt = if isnothing(model.sea_ice) 
        Inf 
    else
        wizard(model.sea_ice)
    end

    atmosphere_Δt = if isnothing(model.atmosphere) 
        Inf 
    else
        wizard(model.atmosphere)
    end

    Δt = min(ocean_Δt, sea_ice_Δt, atmosphere_Δt)

    simulation.Δt = Δt
    
    return nothing
end

