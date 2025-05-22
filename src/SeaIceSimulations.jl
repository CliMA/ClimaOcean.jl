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

function sea_ice_simulation(grid, ocean=nothing;
                            Δt = 5minutes,
                            ice_salinity = 4, # psu
                            advection = nothing, # for the moment
                            tracers = (),
                            ice_heat_capacity = 2100, # J kg⁻¹ K⁻¹
                            ice_consolidation_thickness = 0.05, # m
                            ice_density = 900, # kg m⁻³
                            dynamics = sea_ice_dynamics(grid, ocean),
                            bottom_heat_boundary_condition = nothing,
                            phase_transitions = PhaseTransitions(; ice_heat_capacity, ice_density),
                            conductivity = 2, # kg m s⁻³ K⁻¹
                            internal_heat_flux = ConductiveFlux(; conductivity))

    # Build consistent boundary conditions for the ice model:
    # - bottom -> flux boundary condition
    # - top -> prescribed temperature boundary condition (calculated in the flux computation)

    top_surface_temperature = Field{Center, Center, Nothing}(grid)
    top_heat_boundary_condition = PrescribedTemperature(top_surface_temperature)

    if isnothing(bottom_heat_boundary_condition)
        if isnothing(ocean)
            surface_ocean_salinity = 0
        else
            surface_ocean_salinity = view(ocean.model.tracers.S, :, :, size(ocean.model.grid, 3))
        end
        bottom_heat_boundary_condition = IceWaterThermalEquilibrium(surface_ocean_salinity)
    end

    ice_thermodynamics = SlabSeaIceThermodynamics(grid;
                                                  internal_heat_flux,
                                                  phase_transitions,
                                                  top_heat_boundary_condition,
                                                  bottom_heat_boundary_condition)

    bottom_heat_flux = Field{Center, Center, Nothing}(grid)
    top_heat_flux    = Field{Center, Center, Nothing}(grid)

    # Build the sea ice model
    sea_ice_model = SeaIceModel(grid;
                                ice_salinity,
                                advection,
                                tracers,
                                ice_consolidation_thickness,
                                ice_thermodynamics,
                                dynamics,
                                bottom_heat_flux,
                                top_heat_flux)

    verbose = false

    # Build the simulation
    sea_ice = Simulation(sea_ice_model; Δt, verbose)

    return sea_ice
end

function sea_ice_dynamics(grid, ocean=nothing;
                          sea_ice_ocean_drag_coefficient = 5.5e-3,
                          rheology = ElastoViscoPlasticRheology(pressure_formulation = IceStrength()),
                          coriolis = nothing,
                          solver = SplitExplicitSolver(120))

    if isnothing(ocean)
        SSU = Oceananigans.Fields.ZeroField()
        SSV = Oceananigans.Fields.ZeroField()
    else
        SSU = view(ocean.model.velocities.u, :, :, grid.Nz)
        SSV = view(ocean.model.velocities.v, :, :, grid.Nz)
        if isnothing(coriolis)
            coriolis = ocean.model.coriolis
        end
    end

    τo  = SemiImplicitStress(uₑ=SSU, vₑ=SSV, Cᴰ=sea_ice_ocean_drag_coefficient)
    τua = Field{Face, Center, Nothing}(grid)
    τva = Field{Center, Face, Nothing}(grid)

    return SeaIceMomentumEquation(grid;
                                  coriolis,
                                  top_momentum_stress = (u=τua, v=τva),
                                  bottom_momentum_stress = τo,
                                  rheology,
                                  solver)
end

end
