using Oceananigans
using Oceananigans.Models: update_model_field_time_series!
using Oceananigans.TimeSteppers: Clock
using Oceananigans: SeawaterBuoyancy

using SeawaterPolynomials: TEOS10EquationOfState

import Thermodynamics as AtmosphericThermodynamics  

using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: ClasiusClapyeronSaturation

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!, time
import Oceananigans.Utils: prettytime
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker

struct OceanSeaIceModel{I, A, O, F, C, G} <: AbstractModel{Nothing}
    clock :: C
    grid :: G # TODO: make it so Oceananigans.Simulation does not require this
    atmosphere :: A
    sea_ice :: I
    ocean :: O
    fluxes :: F
end

const OSIM = OceanSeaIceModel

function Base.summary(model::OSIM)
    A = nameof(typeof(architecture(model.grid)))
    G = nameof(typeof(model.grid))
    return string("OceanSeaIceModel{$A, $G}",
                  "(time = ", prettytime(model.clock.time), ", iteration = ", model.clock.iteration, ")")
end

function Base.show(io::IO, cm::OSIM)
    print(io, summary(cm), "\n")
    print(io, "├── ocean: ", summary(cm.ocean.model), "\n")
    print(io, "├── atmosphere: ", summary(cm.atmosphere), "\n")
    print(io, "└── sea_ice: ", summary(cm.sea_ice), "\n")
    return nothing
end

prettytime(model::OSIM)             = prettytime(model.clock.time)
iteration(model::OSIM)              = model.clock.iteration
timestepper(::OSIM)                 = nothing
reset!(::OSIM)                      = nothing
initialize!(::OSIM)                 = nothing
default_included_properties(::OSIM) = tuple()
prognostic_fields(cm::OSIM)         = nothing
fields(::OSIM)                      = NamedTuple()
default_clock(TT)                   = Oceananigans.TimeSteppers.Clock{TT}(0, 0, 1)

reference_density(unsupported) =
    throw(ArgumentError("Cannot extract reference density from $(typeof(unsupported))"))

heat_capacity(unsupported) =
    throw(ArgumentError("Cannot deduce the heat capacity from $(typeof(unsupported))"))

reference_density(ocean::Simulation) = reference_density(ocean.model.buoyancy.formulation)
reference_density(buoyancy_formulation::SeawaterBuoyancy) = reference_density(buoyancy_formulation.equation_of_state)
reference_density(eos::TEOS10EquationOfState) = eos.reference_density
reference_density(sea_ice::SeaIceSimulation) = sea_ice.model.ice_thermodynamics.phase_transitions.ice_density

heat_capacity(ocean::Simulation) = heat_capacity(ocean.model.buoyancy.formulation)
heat_capacity(buoyancy_formulation::SeawaterBuoyancy) = heat_capacity(buoyancy_formulation.equation_of_state)
heat_capacity(sea_ice::SeaIceSimulation) = sea_ice.model.ice_thermodynamics.phase_transitions.ice_heat_capacity

# Does not really matter if there is no model
reference_density(::Nothing) = 0
heat_capacity(::Nothing) = 0

function heat_capacity(::TEOS10EquationOfState{FT}) where FT
    cₚ⁰ = SeawaterPolynomials.TEOS10.teos10_reference_heat_capacity
    return convert(FT, cₚ⁰)
end

function OceanSeaIceModel(ocean, sea_ice=FreezingLimitedOceanTemperature();
                          atmosphere = nothing,
                          radiation = nothing,
                          turbulent_coefficients = nothing,
                          clock = deepcopy(ocean.model.clock),
                          water_vapor_saturation = ClasiusClapyeronSaturation(),
                          ice_vapor_saturation = ClasiusClapyeronSaturation(),
                          water_mole_fraction = 0.98,
                          thermodynamics_parameters = nothing,
                          ocean_reference_density = reference_density(ocean),
                          ocean_heat_capacity = heat_capacity(ocean),
                          ice_reference_density = reference_density(sea_ice),
                          ice_heat_capacity = heat_capacity(sea_ice))

    # Remove some potentially irksome callbacks from the ocean simulation
    pop!(ocean.callbacks, :stop_time_exceeded, nothing)
    pop!(ocean.callbacks, :stop_iteration_exceeded, nothing)
    pop!(ocean.callbacks, :wall_time_limit_exceeded, nothing)
    pop!(ocean.callbacks, :nan_checker, nothing)

    # In case there was any doubt these are meaningless.
    ocean.stop_time = Inf
    ocean.stop_iteration = Inf
    ocean.wall_time_limit = Inf

    if sea_ice isa SeaIceSimulation
        pop!(sea_ice.callbacks, :stop_time_exceeded, nothing)
        pop!(sea_ice.callbacks, :stop_iteration_exceeded, nothing)
        pop!(sea_ice.callbacks, :wall_time_limit_exceeded, nothing)
        pop!(sea_ice.callbacks, :nan_checker, nothing)

        sea_ice.stop_time = Inf
        sea_ice.stop_iteration = Inf
        sea_ice.wall_time_limit = Inf
    end

    # Contains information about flux contributions: bulk formula, prescribed fluxes, etc.
    fluxes = CrossRealmSurfaceFluxes(ocean, sea_ice;
                                     water_vapor_saturation,
                                     ice_vapor_saturation,
                                     water_mole_fraction,
                                     thermodynamics_parameters,
                                     atmosphere, 
                                     turbulent_coefficients,
                                     ocean_reference_density,
                                     ocean_heat_capacity,
                                     ice_reference_density,
                                     ice_heat_capacity,
                                     radiation)

    ocean_sea_ice_model = OceanSeaIceModel(clock,
                                           ocean.model.grid,
                                           atmosphere,
                                           sea_ice,
                                           ocean,
                                           fluxes)

    update_state!(ocean_sea_ice_model)

    return ocean_sea_ice_model
end

time(coupled_model::OceanSeaIceModel) = coupled_model.clock.time

# Check for NaNs in the first prognostic field (generalizes to prescribed velocitries).
function default_nan_checker(model::OceanSeaIceModel)
    u_ocean = model.ocean.model.velocities.u
    nan_checker = NaNChecker((; u_ocean))
    return nan_checker
end
