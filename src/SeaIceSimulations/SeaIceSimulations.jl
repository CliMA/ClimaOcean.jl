module SeaIceSimulations

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: with_tracers
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.BoundaryConditions: DefaultBoundaryCondition
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node
using Oceananigans.OrthogonalSphericalShellGrids
using Oceananigans.Operators

export sea_ice_simulation

import ClimaOcean.OceanSeaIceModels:
    sea_ice_thickness,
    sea_ice_concentration,
    reference_density,
    heat_capacity, 
    interpolate_sea_ice_state!,
    compute_net_sea_ice_fluxes!,
    compute_sea_ice_ocean_fluxes!

import Oceananigans.TimeSteppers: time_step!

include("freezing_limited_ocean_temperature.jl")
include("sea_ice_simulation.jl")
include("assemble_net_sea_ice_fluxes.jl")

# When using an ClimaSeaIce simulation, we assume that the exchange grid is the sea-ice grid
interpolate_sea_ice_state!(interfaces, ::Simulation{<:SeaIceModel}, coupled_model) = nothing
interpolate_sea_ice_state!(interfaces, ::FreezingLimitedOceanTemperature, coupled_model) = nothing

end
