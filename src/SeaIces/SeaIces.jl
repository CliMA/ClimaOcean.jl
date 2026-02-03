module SeaIces

export sea_ice_simulation, FreezingLimitedOceanTemperature

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils
using Oceananigans.Utils: with_tracers
using Oceananigans.Grids: architecture
using Oceananigans.Fields: ZeroField
using Oceananigans.BoundaryConditions: DefaultBoundaryCondition
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node
using Oceananigans.OrthogonalSphericalShellGrids
using Oceananigans.Operators
using KernelAbstractions: @kernel, @index

import ClimaOcean.OceanSeaIceModels: interpolate_state!,
                                     sea_ice_concentration,
                                     sea_ice_thickness,
                                     reference_density,
                                     heat_capacity,
                                     update_net_fluxes!,
                                     default_sea_ice

import ClimaOcean.OceanSeaIceModels.InterfaceComputations: ComponentExchanger,
                                                           compute_atmosphere_sea_ice_fluxes!,
                                                           compute_sea_ice_ocean_fluxes!,
                                                           net_fluxes,
                                                           ThreeEquationHeatFlux,
                                                           default_ai_temperature

import Oceananigans.TimeSteppers: time_step!

include("freezing_limited_ocean_temperature.jl")
include("sea_ice_simulation.jl")
include("assemble_net_sea_ice_fluxes.jl")

default_sea_ice() = FreezingLimitedOceanTemperature()

# When using an ClimaSeaIce simulation, we assume that the exchange grid is the sea-ice grid
interpolate_state!(exchanger, grid, ::Simulation{<:SeaIceModel},       coupled_model) = nothing
interpolate_state!(exchanger, grid, ::FreezingLimitedOceanTemperature, coupled_model) = nothing

# ComponentExchangers
ComponentExchanger(sea_ice::FreezingLimitedOceanTemperature, grid) = nothing

function ComponentExchanger(sea_ice::Simulation{<:SeaIceModel}, grid) 
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

    return ComponentExchanger((; u, v, h, hc, ℵ), nothing)
end

end
