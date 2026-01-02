using ClimaOcean.Oceans

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: 
    net_fluxes,
    sea_ice_ocean_interface, 
    atmosphere_ocean_interface, 
    initialize!,
    ComponentExchanger,
    default_exchange_grid

import ClimaOcean.OceanSeaIceModels:
    interpolate_state!,
    update_net_fluxes!

import ClimaOcean.Oceans: get_radiative_forcing

function ComponentExchanger(ocean::VerosOceanSimulation, grid) 
    state = (; u = Field{Face, Center, Nothing}(grid),
               v = Field{Center, Face, Nothing}(grid),
               T = Field{Center, Center, Nothing}(grid),
               S = Field{Center, Center, Nothing}(grid))

    return ComponentExchanger(state, nothing)
end

default_exchange_grid(atmosphere, ocean::VerosOceanSimulation, sea_ice) = surface_grid(ocean)

@inline function net_fluxes(ocean::VerosOceanSimulation)
    grid = surface_grid(ocean)
    u = Field{Face,   Center, Nothing}(grid)
    v = Field{Center, Face,   Nothing}(grid)
    T = Field{Center, Center, Nothing}(grid)
    S = Field{Center, Center, Nothing}(grid)

    return (; u, v, T, S)
end

function interpolate_state!(exchanger, exchange_grid, ocean::VerosOceanSimulation, coupled_model)
    u = exchanger.state.u
    v = exchanger.state.v
    T = exchanger.state.T
    S = exchanger.state.S

    set!(u, ocean.setup.state.variables.u)
    set!(v, ocean.setup.state.variables.v)
    set!(T, ocean.setup.state.variables.temp)
    set!(S, ocean.setup.state.variables.salt)

    return nothing
end

initialize!(exchanger::ComponentExchanger, grid, ::VerosOceanSimulation) = nothing

get_radiative_forcing(ocean::VerosOceanSimulation) = nothing

function update_net_fluxes!(coupled_model, ocean::VerosOceanSimulation)

    # Update the flux containers
    Oceans.update_net_ocean_fluxes!(coupled_model, ocean, coupled_model.interfaces.exchanger.grid)
    net_ocean_fluxes = coupled_model.interfaces.net_fluxes.ocean
   
    # Pass the flux values to the python ocean
    nx = pyconvert(Int, ocean.setup.state.settings.nx) + 4
    ny = pyconvert(Int, ocean.setup.state.settings.ny) + 4

    ρₒ = pyconvert(eltype(ocean), ocean.setup.state.settings.rho_0)
    taux = view(parent(net_ocean_fluxes.u), 1:nx, 1:ny, 1) .* ρₒ
    tauy = view(parent(net_ocean_fluxes.v), 1:nx, 1:ny, 1) .* ρₒ

    set!(ocean, "surface_taux", taux; path=:variables)
    set!(ocean, "surface_tauy", tauy; path=:variables)

    temp_flux = view(parent(net_ocean_fluxes.T), 1:nx, 1:ny, 1)
    salt_flux = view(parent(net_ocean_fluxes.S), 1:nx, 1:ny, 1)

    set!(ocean, "forc_temp_surface", temp_flux; path=:variables)
    set!(ocean, "forc_salt_surface", salt_flux; path=:variables)

    return nothing
end
