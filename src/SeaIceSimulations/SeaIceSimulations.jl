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
    compute_net_sea_ice_fluxes!,
    compute_sea_ice_ocean_fluxes!

import Oceananigans.TimeSteppers: time_step!

include("freezing_limited_ocean_temperature.jl")
include("sea_ice_simulation.jl")

end
