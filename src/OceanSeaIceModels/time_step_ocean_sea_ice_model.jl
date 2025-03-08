using .InterfaceComputations:
    compute_atmosphere_ocean_fluxes!,
    compute_sea_ice_ocean_fluxes!,
    compute_net_ocean_fluxes!,
    compute_net_sea_ice_fluxes!,
    interpolate_atmosphere_state!

using ClimaSeaIce: SeaIceModel, SeaIceThermodynamics
using Oceananigans.Grids: φnode

using Printf

function time_step!(coupled_model::OceanSeaIceModel, Δt; callbacks=[], compute_tendencies=true)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    atmosphere = coupled_model.atmosphere
    clock = coupled_model.clock

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

        time_step!(sea_ice, Δt)
    end

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
