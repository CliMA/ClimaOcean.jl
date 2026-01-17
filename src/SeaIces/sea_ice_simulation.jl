using ClimaSeaIce
using ClimaSeaIce: SeaIceModel, SlabSeaIceThermodynamics, PhaseTransitions, ConductiveFlux
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver, SemiImplicitStress, SeaIceMomentumEquation, StressBalanceFreeDrift
using ClimaSeaIce.Rheologies: IceStrength, ElastoViscoPlasticRheology

using ClimaOcean.OceanSeaIceModels: ocean_surface_salinity, ocean_surface_velocities
using ClimaOcean.Oceans: Default

default_rotation_rate = Oceananigans.defaults.planet_rotation_rate

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
                            top_heat_boundary_condition = nothing,
                            phase_transitions = PhaseTransitions(; ice_heat_capacity, ice_density),
                            conductivity = 2, # kg m s⁻³ K⁻¹
                            internal_heat_flux = ConductiveFlux(; conductivity))

    # Build consistent boundary conditions for the ice model:
    # - bottom -> flux boundary condition
    # - top -> prescribed temperature boundary condition (calculated in the flux computation)

    if isnothing(top_heat_boundary_condition)
        top_surface_temperature = Field{Center, Center, Nothing}(grid)
        top_heat_boundary_condition = PrescribedTemperature(top_surface_temperature.data)
    end

    if isnothing(bottom_heat_boundary_condition)
        if isnothing(ocean)
            surface_ocean_salinity = 0
        else
            kᴺ = size(grid, 3)
            surface_ocean_salinity = ocean_surface_salinity(ocean)
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
                          rheology = ElastoViscoPlasticRheology(),
                          coriolis = HydrostaticSphericalCoriolis(; rotation_rate=default_rotation_rate),
                          free_drift = nothing,
                          solver = SplitExplicitSolver(120))

    SSU, SSV = ocean_surface_velocities(ocean)
    sea_ice_ocean_drag_coefficient = convert(eltype(grid), sea_ice_ocean_drag_coefficient)

    τo  = SemiImplicitStress(uₑ=SSU, vₑ=SSV, Cᴰ=sea_ice_ocean_drag_coefficient)
    τua = Field{Face, Center, Nothing}(grid)
    τva = Field{Center, Face, Nothing}(grid)

    if isnothing(free_drift)
        free_drift = StressBalanceFreeDrift((u=τua, v=τva), τo)
    end

    return SeaIceMomentumEquation(grid;
                                  coriolis,
                                  top_momentum_stress = (u=τua, v=τva),
                                  bottom_momentum_stress = τo,
                                  rheology,
                                  free_drift,
                                  solver)
end

#####
##### Extending OceanSeaIceModels interface
#####

sea_ice_thickness(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thickness
sea_ice_concentration(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_concentration

heat_capacity(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thermodynamics.phase_transitions.ice_heat_capacity
reference_density(sea_ice::Simulation{<:SeaIceModel}) = sea_ice.model.ice_thermodynamics.phase_transitions.ice_density

function net_fluxes(sea_ice::Simulation{<:SeaIceModel}) 
    net_momentum_fluxes = if isnothing(sea_ice.model.dynamics)
        u = Field{Face, Center, Nothing}(sea_ice.model.grid)
        v = Field{Center, Face, Nothing}(sea_ice.model.grid)
        (; u, v)
    else
        u = sea_ice.model.dynamics.external_momentum_stresses.top.u
        v = sea_ice.model.dynamics.external_momentum_stresses.top.v
        (; u, v)
    end

    net_top_sea_ice_fluxes = merge((; heat=sea_ice.model.external_heat_fluxes.top), net_momentum_fluxes)
    net_bottom_sea_ice_fluxes = (; heat=sea_ice.model.external_heat_fluxes.bottom)

    return (; bottom = net_bottom_sea_ice_fluxes, top = net_top_sea_ice_fluxes)
end

function default_ai_temperature(sea_ice::Simulation{<:SeaIceModel})
    conductive_flux = sea_ice.model.ice_thermodynamics.internal_heat_flux.parameters.flux
    return SkinTemperature(conductive_flux)
end

# Constructor that accepts the sea-ice model
function ThreeEquationHeatFlux(sea_ice::Simulation{<:SeaIceModel}, FT::DataType = Oceananigans.defaults.FloatType; 
                               heat_transfer_coefficient = 0.0095,
                               salt_transfer_coefficient = heat_transfer_coefficient / 35,
                               friction_velocity = convert(FT, 0.002))

    conductive_flux = sea_ice.model.ice_thermodynamics.internal_heat_flux.parameters.flux
    ice_temperature = sea_ice.model.ice_thermodynamics.top_surface_temperature
    
    return ThreeEquationHeatFlux(conductive_flux,
                                 ice_temperature,
                                 convert(FT, heat_transfer_coefficient),
                                 convert(FT, salt_transfer_coefficient),
                                 friction_velocity) 
end
