module SeaIceSimulations

export sea_ice_simulation

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: with_tracers
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.BoundaryConditions: DefaultBoundaryCondition
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node
using Oceananigans.OrthogonalSphericalShellGrids
using Oceananigans.Operators

using ClimaOcean.OceanSeaIceModels: OceanSeaIceModel
using ClimaOcean: OceanSeaIceModels

import Oceananigans.TimeSteppers: time_step!

include("freezing_limited_ocean_temperature.jl")
include("sea_ice_simulation.jl")
include("assemble_net_sea_ice_fluxes.jl")

const NoSeaIceModel = Union{OceanSeaIceModel{Nothing}, OceanSeaIceModel{<:FreezingLimitedOceanTemperature}}

# When using an ClimaSeaIce simulation, we assume that the exchange grid is the sea-ice grid
OceanSeaIceModels.interpolate_sea_ice_state!(interfaces, ::Simulation{<:SeaIceModel}, coupled_model) = nothing
OceanSeaIceModels.interpolate_sea_ice_state!(interfaces, ::FreezingLimitedOceanTemperature, coupled_model) = nothing
OceanSeaIceModels.ComponentExchanger(sea_ice::FreezingLimitedOceanTemperature, grid) = nothing

function OceanSeaIceModels.ComponentExchanger(sea_ice::Simulation{<:SeaIceModel}, grid) 
    sea_ice_grid = sea_ice.model.grid
    
    if sea_ice_grid == grid
        u  = sea_ice.model.velocities.u
        v  = sea_ice.model.velocities.v
        h  = sea_ice.model.ice_thickness
        hc = sea_ice.model.ice_consolidation_thickness
        ℵ  = sea_ice.model.ice_concentration
    else
        u  = Field{Center, Center, Nothing}(grid)
        v  = Field{Center, Center, Nothing}(grid)
        h  = Field{Center, Center, Nothing}(grid)
        hc = Field{Center, Center, Nothing}(grid)
        ℵ  = Field{Center, Center, Nothing}(grid)
    end

    return OceanSeaIceModels.ComponentExchanger((; u, v, h, hc, ℵ), nothing)
end

end
