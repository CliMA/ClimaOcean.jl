using .InterfaceComputations:
    compute_atmosphere_ocean_fluxes!,
    compute_sea_ice_ocean_fluxes!,
    compute_net_ocean_fluxes!,
    compute_net_sea_ice_fluxes!,
    interpolate_atmosphere_state!

using ClimaSeaIce: SeaIceModel, SeaIceThermodynamics
using Oceananigans.Grids: φnode
using Oceananigans.Simulations: Callback, TimeStepCallsite

using Printf

#####
##### Component time step utilities
#####

# Get time step from component simulation, or return Inf if not applicable
component_Δt(component) = Inf
component_Δt(sim::Simulation) = sim.Δt
component_Δt(::PrescribedAtmosphere) = Inf  # prescribed atmosphere has no CFL constraint
component_Δt(::FreezingLimitedOceanTemperature) = Inf  # not a dynamical model

"""
    align_component_steps!(simulation)

Synchronize the coupled simulation's time step with its component simulations.

This function:
1. Invokes callbacks on component simulations (e.g., TimeStepWizard updates ocean.Δt)
2. Sets the coupled simulation's Δt to the minimum of all component time steps

This is useful when component simulations have adaptive time-stepping (via TimeStepWizard)
and you want the coupled model to automatically use the most restrictive time step.

Example
=======

```julia
ocean = nonhydrostatic_ocean_simulation(grid)  # has TimeStepWizard
coupled_model = OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model; Δt=10.0, stop_time=1hour)

# Add callback to synchronize time steps before each coupled time step
add_callback!(simulation, align_component_steps!)
```
"""
function align_component_steps!(simulation)
    coupled_model = simulation.model
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    atmosphere = coupled_model.atmosphere

    # Set coupled simulation Δt to minimum of component time steps.
    # Note: we don't include the current simulation.Δt in the min, so the
    # time step can grow as well as shrink based on component constraints.
    synced_Δt = min(component_Δt(ocean),
                    component_Δt(sea_ice),
                    component_Δt(atmosphere))

    if isfinite(synced_Δt)
        simulation.Δt = synced_Δt
    end

    return nothing
end

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
