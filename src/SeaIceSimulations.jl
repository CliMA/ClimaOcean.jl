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
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies

using ClimaOcean.OceanSimulations: Default

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

function default_sea_ice_dynamics(grid; 
                                  ocean, # Cannot do it without an ocean
                                  sea_ice_ocean_drag_coefficient = 5.5e-3,
                                  rheology = ElastoViscoPlasticRheology(),
                                  solver = SplitExplicitSolver(120))

    SSU = view(ocean.model.velocities.u, :, :, grid.Nz)
    SSV = view(ocean.model.velocities.v, :, :, grid.Nz)

    τo  = SemiImplicitStress(uₑ=SSU, vₑ=SSV, Cᴰ=sea_ice_ocean_drag_coefficient)
    τua = Field{Face, Center, Nothing}(grid)
    τva = Field{Center, Face, Nothing}(grid)

    return SeaIceMomentumEquation(grid;
                                  coriolis = ocean.model.coriolis,
                                  top_momentum_stress = (u=τua, v=τva),
                                  bottom_momentum_stress = τo,
                                  ocean_velocities = (u=0.1*SSU, v=0.1*SSV),
                                  rheology,
                                  solver)
end

end