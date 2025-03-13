module SeaIceSimulations

export sea_ice_simulation

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: with_tracers
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.BoundaryConditions: DefaultBoundaryCondition
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node
using Oceananigans.OrthogonalSphericalShellGrids

using Oceananigans.BuoyancyFormulations: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Operators

using ClimaSeaIce
using ClimaSeaIce: SeaIceModel, SlabSeaIceThermodynamics, PhaseTransitions, ConductiveFlux
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium

using ClimaOcean.OceanSimulations: Default

import ClimaOcean: reference_density, heat_capacity

reference_density(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thermodynamics.phase_transitions.ice_density
heat_capacity(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thermodynamics.phase_transitions.ice_heat_capacity

function sea_ice_simulation(grid;
                            Δt = 5minutes,
                            ice_salinity = 0, # psu
                            advection = nothing, # for the moment
                            tracers = (), 
                            ice_heat_capacity = 2100, # J kg⁻¹ K⁻¹
                            ice_consolidation_thickness = 0.05, # m
                            ice_density = 900, # kg m⁻³
                            dynamics = nothing,
                            bottom_heat_boundary_condition = IceWaterThermalEquilibrium(),
                            phase_transitions = PhaseTransitions(; ice_heat_capacity, ice_density),
                            conductivity = 2, # kg m s⁻³ K⁻¹
                            internal_heat_flux = ConductiveFlux(; conductivity))

    # Build consistent boundary conditions for the ice model:
    # - bottom -> flux boundary condition
    # - top -> prescribed temperature boundary condition (calculated in the flux computation)

    top_surface_temperature = Field{Center, Center, Nothing}(grid)
    top_heat_boundary_condition = PrescribedTemperature(top_surface_temperature)

    ice_thermodynamics = SlabSeaIceThermodynamics(grid;
                                                  internal_heat_flux,
                                                  phase_transitions,
                                                  top_heat_boundary_condition,
                                                  bottom_heat_boundary_condition)

    bottom_heat_flux = Field{Center, Center, Nothing}(grid)
    top_heat_flux    = Field{Center, Center, Nothing}(grid)

    # top_momentum_stress = (u = Field{Face, Center, Nothing}(grid),
    #                        v = Field{Center, Face, Nothing}(grid))

    # Build the sea ice model
    sea_ice_model = SeaIceModel(grid;
                                ice_salinity,
                                advection,
                                tracers,
                                ice_consolidation_thickness,
                                # top_momentum_stress,
                                ice_thermodynamics,
                                dynamics,
                                bottom_heat_flux,
                                top_heat_flux)

    verbose = false

    # Build the simulation
    sea_ice = Simulation(sea_ice_model; Δt, verbose)

    return sea_ice
end

end