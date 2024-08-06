module SeaIceSimulations

export sea_ice_simulation

using Oceananigans
using Oceananigans.Units
using Oceananigans.Advection: TracerAdvection
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node

using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: SlabSeaIceThermodynamics
using ClimaSeaIce.SeaIceDynamics: ExplicitMomentumSolver

using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Operators

default_sea_ice_advection() = WENO(; order = 7)

function sea_ice_simulation(grid; Δt = 5minutes,
                            reference_density = 1020,
                            rotation_rate = Ω_Earth,
                            gravitational_acceleration = g_Earth,
                            ocean_ice_drag_coefficient = 0.0055,
                            ocean_velocities = nothing,
                            top_heat_flux = nothing,
                            ice_salinity = 4,
                            ice_thermodynamics = SlabSeaIceThermodynamics(grid),
                            ice_dynamics = ExplicitMomentumSolver(grid; substeps = 100, ocean_ice_drag_coefficient),
                            coriolis = HydrostaticSphericalCoriolis(; rotation_rate),
                            advection = default_sea_ice_advection(),
                            verbose = false)

    # Set up boundary conditions using Field
    top_zonal_momentum_flux      = Jᵘ = Field{Face, Center, Nothing}(grid)
    top_meridional_momentum_flux = Jᵛ = Field{Center, Face, Nothing}(grid)

    Jᵀ = isnothing(top_heat_flux) ? Field{Center, Center, Nothing}(grid) : top_heat_flux

    sea_ice_model = SeaIceModel(grid;
                                advection,
                                ice_thermodynamics, 
                                ice_dynamics, 
                                ocean_velocities,
                                top_u_stress = Jᵘ,
                                top_v_stress = Jᵛ,
                                coriolis,
                                top_heat_flux = Jᵀ,
                                ice_salinity)

    sea_ice = Simulation(sea_ice_model; Δt, verbose)

    return sea_ice
end

end # module
