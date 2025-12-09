
struct StateExchanger{G, AS, OS, SS, AEX, OEX, SIE}
    grid :: G
    atmosphere_exchanger :: AST
    ocean_exchanger :: OS
    sea_ice_exchanger :: SS
end

struct AtmosphereExchanger{FC, CF, CC, AEX}
    u  :: FC
    v  :: CF
    T  :: CC
    q  :: CC
    p  :: CC
    Qs :: CC
    Qℓ :: CC
    Mp :: CC
    exchanger :: AEX
end

struct OceanExchanger{FC, CF, CC, OEX}
    u :: FC
    v :: CF
    T :: CC
    S :: CC
    exchanger :: OEX
end

struct SeaIceExchanger{FC, CF, CC, SIE}
    u :: FC
    v :: CF
    h :: CC
    ℵ :: CC
    exchanger :: SIE
end

function StateExchanger(grid, ocean, atmosphere, sea_ice)
    # TODO: generalize this
    atmosphere_state = ExchangeAtmosphereState(atmosphere, grid)
    ocean_state = ExchangeOceanState(ocean, grid)
    sea_ice_state = ExchangeSeaIceState(sea_ice, grid)

    return StateExchanger(grid, atmosphere_state, ocean_state, sea_ice_state)
end

function initialize!(exchanger::StateExchanger, model)
    initialize!(exchanger.atmosphere_exchanger, model.atmosphere)
    initialize!(exchanger.ocean_exchanger,      model.ocean)
    initialize!(exchanger.sea_ice_exchanger,    model.sea_ice)
    return nothing
end
