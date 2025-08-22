using Oceananigans.Models: initialization_update_state!

using ClimaOcean.OceanSeaIceModels.InterfaceComputations: ExchangeAtmosphereState, 
                                                          atmosphere_exchanger, 
                                                          SimilarityTheoryFluxes, 
                                                          Radiation

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: 
    state_exchanger, 
    sea_ice_ocean_interface, 
    atmosphere_ocean_interface, 
    initialize!,
    get_ocean_state,
    ocean_surface_fluxes,
    get_radiative_forcing,
    fill_up_net_fluxes!

mutable struct VerosStateExchanger{G, OST, AST, AEX}
    exchange_grid :: G
    exchange_ocean_state :: OST
    exchange_atmosphere_state :: AST
    atmosphere_exchanger :: AEX
end

mutable struct ExchangeOceanState{FC, CF, CC}
    u  :: FC
    v  :: CF
    T  :: CC
    S  :: CC
end

ExchangeOceanState(grid) = ExchangeOceanState(Field{Face, Center, Nothing}(grid),
                                              Field{Center, Face, Nothing}(grid),
                                              Field{Center, Center, Nothing}(grid),
                                              Field{Center, Center, Nothing}(grid))

function state_exchanger(ocean::VerosOceanSimulation, atmosphere)
    exchange_grid = surface_grid(ocean)
    exchange_ocean_state = ExchangeOceanState(exchange_grid)
    exchange_atmosphere_state = ExchangeAtmosphereState(exchange_grid)

    atmos_exchanger = atmosphere_exchanger(atmosphere, exchange_grid)

    return VerosStateExchanger(exchange_grid,
                               exchange_ocean_state,
                               exchange_atmosphere_state,
                               atmos_exchanger)
end

atmosphere_ocean_interface(ocean::VerosOceanSimulation, args...) = 
    atmosphere_ocean_interface(surface_grid(ocean), args...)

sea_ice_ocean_interface(ocean::VerosOceanSimulation, args...) = 
    sea_ice_ocean_interface(surface_grid(ocean), args...)

initialize!(exchanger::VerosStateExchanger, atmosphere) = 
    initialize!(exchanger.atmosphere_exchanger, exchanger.exchange_grid, atmosphere)

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

@inline function ocean_surface_fluxes(ocean::VerosOceanSimulation)
    grid = surface_grid(ocean)
    u = Field{Face, Center, Nothing}(grid)
    v = Field{Center, Face, Nothing}(grid)
    T = Field{Center, Center, Nothing}(grid)
    S = Field{Center, Center, Nothing}(grid)

    return (; u, v, T, S)
end

@inline get_radiative_forcing(ocean::VerosOceanSimulation) = nothing

function fill_net_fluxes!(ocean::VerosOceanSimulation, net_ocean_fluxes)
    nx = pyconvert(Int, ocean.setup.state.settings.nx) + 4
    ny = pyconvert(Int, ocean.setup.state.settings.ny) + 4
    t1 = view(parent(net_ocean_fluxes.u), 1:nx, 1:ny, 1)
    t2 = view(parent(net_ocean_fluxes.v), 1:nx, 1:ny, 1)

    # TODO: Remove this when they fix veros to
    # be able to force with an array instead of Climatology
    ta = zeros(size(t1)..., 12)
    tb = zeros(size(t2)..., 12)
    for t in 1:12
        ta[:, :, t] .= t1
        tb[:, :, t] .= t2
    end

    veros_set!(ocean, "taux", ta)
    veros_set!(ocean, "tauy", tb)

    # TODO: Add heat flux and salinity flux when they
    # fix veros to be able to force with prescribed boundary
    # conditions rather than restoring

    # t1 = view(parent(net_ocean_fluxes.T), 1:nx, 2:ny+1, 1)
    # t2 = view(parent(net_ocean_fluxes.S), 1:nx, 2:ny+1, 1)

    # veros_set!(ocean, "temp_flux", t1)
    # veros_set!(ocean, "salt_flux", t2)

    return nothing
end
