struct StateExchanger{G, AS, OS, SS, AEX, OEX, SIE}
    grid :: G
    atmosphere :: AST
    ocean :: OS
    sea_ice :: SS
end

struct ComponentExchanger{AS, AEX}
    state :: AS
    exchanger :: AEX
end

# For ``nothing'' components, we don't need an exchanger
ComponentExchanger(::Nothing, grid) = nothing

function StateExchanger(grid, ocean, atmosphere, sea_ice)
    # TODO: generalize this
    atmosphere_state = ComponentExchanger(atmosphere, grid)
    ocean_state      = ComponentExchanger(ocean, grid)
    sea_ice_state    = ComponentExchanger(sea_ice, grid)

    return StateExchanger(grid, atmosphere_state, ocean_state, sea_ice_state)
end

function initialize!(exchanger::StateExchanger, model)
    initialize!(exchanger.atmosphere, exchanger.grid, model.atmosphere)
    initialize!(exchanger.ocean,      exchanger.grid, model.ocean)
    initialize!(exchanger.sea_ice,    exchanger.grid, model.sea_ice)
    return nothing
end
