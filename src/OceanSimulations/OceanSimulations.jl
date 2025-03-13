module OceanSimulations

export ocean_simulation, PrescribedOcean, ocean_state

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: with_tracers
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.BoundaryConditions: DefaultBoundaryCondition
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node, MutableGridOfSomeKind
using Oceananigans.OrthogonalSphericalShellGrids

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    CATKEMixingLength,
    CATKEEquation

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using Oceananigans.BuoyancyFormulations: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Operators

import ClimaOcean: reference_density, heat_capacity

reference_density(ocean::Simulation{<:HydrostaticFreeSurfaceModel}) = reference_density(ocean.model.buoyancy.formulation)
reference_density(buoyancy_formulation::SeawaterBuoyancy) = reference_density(buoyancy_formulation.equation_of_state)
reference_density(eos::TEOS10EquationOfState) = eos.reference_density

heat_capacity(ocean::Simulation{<:HydrostaticFreeSurfaceModel}) = heat_capacity(ocean.model.buoyancy.formulation)
heat_capacity(buoyancy_formulation::SeawaterBuoyancy) = heat_capacity(buoyancy_formulation.equation_of_state)

struct Default{V}
    value :: V
end

"""
    default_or_override(default::Default, alternative_default=default.value) = alternative_default
    default_or_override(override, alternative_default) = override

Either return `default.value`, an `alternative_default`, or an `override`.

The purpose of this function is to help define constructors with "configuration-dependent" defaults.
For example, the default bottom drag should be 0 for a single column model, but 0.003 for a global model.
We therefore need a way to specify both the "normal" default 0.003 as well as the "alternative default" 0,
all while respecting user input and changing this to a new value if specified.
"""
default_or_override(default::Default, possibly_alternative_default=default.value) = possibly_alternative_default
default_or_override(override, alternative_default=nothing) = override

include("barotropic_potential_forcing.jl")
include("ocean_simulation.jl")
include("prescribed_ocean.jl")

# Grab the ocean state (Needs extending for new ocean models)
ocean_state(::Nothing) = (u = nothing, v = nothing, T = nothing, S = nothing)
ocean_state(ocean::Simulation{<:HydrostaticFreeSurfaceModel}) = ocean_state(ocean.model)
ocean_state(ocean) = (u = ocean.velocities.u,
                      v = ocean.velocities.v,
                      T = ocean.tracers.T,
                      S = ocean.tracers.S)

end # module
