using Oceananigans.Models: initialization_update_state!

using ClimaOcean.OceanSeaIceModels.InterfaceComputations: ExchangeAtmosphereState, 
                                                          atmosphere_exchanger, 
                                                          SimilarityTheoryFluxes, 
                                                          Radiation

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: 
    sea_ice_ocean_interface, 
    atmosphere_ocean_interface, 
    initialize!,
    ComponentExchanger,
    default_exchange_grid

import ClimaOcean.Oceans:
    ocean_salinity,
    ocean_temperature,
    ocean_suface_salinity,
    ocean_surface_temperature,
    ocean_suface_velocity,
    get_radiative_forcing

function ComponentExchanger(ocean::VerosOceanSimulation, grid) 
    state = (; u  = Field{Face, Center, Nothing}(grid),
               v  = Field{Center, Face, Nothing}(grid),
               T  = Field{Center, Center, Nothing}(grid),
               S  = Field{Center, Center, Nothing}(grid))

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

initialize!(exchanger, grid, ::VerosOceanSimulation) = nothing

@inline function get_ocean_state(ocean::VerosOceanSimulation, exchanger::VerosStateExchanger)
    u = exchanger.exchange_ocean_state.u
    v = exchanger.exchange_ocean_state.v
    T = exchanger.exchange_ocean_state.T
    S = exchanger.exchange_ocean_state.S

    set!(u, ocean.setup.state.variables.u)
    set!(v, ocean.setup.state.variables.v)
    set!(T, ocean.setup.state.variables.temp)
    set!(S, ocean.setup.state.variables.salt)

    return (; u, v, T, S)
end

@inline get_radiative_forcing(ocean::VerosOceanSimulation) = nothing

function fill_net_fluxes!(ocean::VerosOceanSimulation, net_ocean_fluxes)
    nx = pyconvert(Int, ocean.setup.state.settings.nx) + 4
    ny = pyconvert(Int, ocean.setup.state.settings.ny) + 4

    ρₒ = pyconvert(eltype(ocean), ocean.setup.state.settings.rho_0)
    taux = view(parent(net_ocean_fluxes.u), 1:nx, 1:ny, 1) .* ρₒ
    tauy = view(parent(net_ocean_fluxes.v), 1:nx, 1:ny, 1) .* ρₒ

    set!(ocean, "surface_taux", tx; path=:variables)
    set!(ocean, "surface_tauy", ty; path=:variables)

    temp_flux = view(parent(net_ocean_fluxes.T), 1:nx, 1:ny, 1)
    salt_flux = view(parent(net_ocean_fluxes.S), 1:nx, 1:ny, 1)

    set!(ocean, "forc_temp_surface", temp_flux; path=:variables)
    set!(ocean, "forc_salt_surface", salt_flux; path=:variables)

    return nothing
end
