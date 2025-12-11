struct ComponentExchanger{S, EX}
    state :: S
    regridder :: EX
end

struct StateExchanger{G, A, O, S}
    grid :: G
    atmosphere :: A
    ocean :: O
    sea_ice :: S

    function StateExchanger(grid, atmosphere, ocean, sea_ice)
        atmosphere_exchanger = ComponentExchanger(atmosphere, grid)
        ocean_exchanger      = ComponentExchanger(ocean, grid)
        sea_ice_exchanger    = ComponentExchanger(sea_ice, grid)

        G = typeof(grid)
        A = typeof(atmosphere_exchanger)
        O = typeof(ocean_exchanger)
        S = typeof(sea_ice_exchanger)
        
        return new{G, A, O, S}(grid, 
                               atmosphere_exchanger, 
                               ocean_exchanger, 
                               sea_ice_exchanger)
    end
end

# For ``nothing'' components, we don't need an exchanger
ComponentExchanger(::Nothing, grid) = nothing

function initialize!(exchanger::StateExchanger, model)
    initialize!(exchanger.atmosphere, exchanger.grid, model.atmosphere)
    initialize!(exchanger.ocean,      exchanger.grid, model.ocean)
    initialize!(exchanger.sea_ice,    exchanger.grid, model.sea_ice)
    return nothing
end

# fallback
initialize!(::Nothing, grid, component) = nothing
initialize!(exchanger::ComponentExchanger, grid, component) = nothing
