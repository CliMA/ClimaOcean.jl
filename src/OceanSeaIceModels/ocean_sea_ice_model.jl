using Oceananigans
using Oceananigans.TimeSteppers: Clock
using Oceananigans: SeawaterBuoyancy
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using KernelAbstractions: @kernel, @index
using SeawaterPolynomials: TEOS10EquationOfState

import Thermodynamics as AtmosphericThermodynamics

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Architectures: architecture
import Oceananigans.Fields: set!
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker, initialization_update_state!
import Oceananigans.OutputWriters: default_included_properties, checkpointer_address,
                                   write_output!, initialize_jld2_file!,
                                   required_checkpointed_properties, default_checkpointed_properties
import Oceananigans.Simulations: reset!, initialize!, iteration, run!
import Oceananigans.TimeSteppers: time_step!, update_state!, time
import Oceananigans.Utils: prettytime

import .PrescribedAtmospheres: set_clock!

mutable struct OceanSeaIceModel{I, A, O, F, C, Arch} <: AbstractModel{Nothing, Arch}
    architecture :: Arch
    clock :: C
    atmosphere :: A
    sea_ice :: I
    ocean :: O
    interfaces :: F
end

const OSIM = OceanSeaIceModel
const OSIMPA = OceanSeaIceModel{<:Any, <:PrescribedAtmosphere}

function Base.summary(model::OSIM)
    A = nameof(typeof(architecture(model)))
    return string("OceanSeaIceModel{$A}",
                  "(time = ", prettytime(model.clock.time), ", iteration = ", model.clock.iteration, ")")
end

function Base.show(io::IO, cm::OSIM)

    if cm.sea_ice isa Simulation
        sea_ice_summary = summary(cm.sea_ice.model)
    else
        sea_ice_summary = summary(cm.sea_ice)
    end

    print(io, summary(cm), "\n")
    print(io, "├── ocean: ", summary(cm.ocean.model), "\n")
    print(io, "├── atmosphere: ", summary(cm.atmosphere), "\n")
    print(io, "├── sea_ice: ", sea_ice_summary, "\n")
    print(io, "└── interface: ", summary(cm.interfaces))
    return nothing
end

# Assumption: We have an ocean!
architecture(model::OSIM)                = architecture(model.ocean.model)
Base.eltype(model::OSIM)                 = Base.eltype(model.ocean.model)
prettytime(model::OSIM)                  = prettytime(model.clock.time)
iteration(model::OSIM)                   = model.clock.iteration
timestepper(::OSIM)                      = nothing
default_included_properties(::OSIM)      = tuple()
prognostic_fields(::OSIM)                = nothing
fields(::OSIM)                           = NamedTuple()
default_clock(TT)                        = Oceananigans.TimeSteppers.Clock{TT}(0, 0, 1)
time(model::OSIM)                        = model.clock.time
required_checkpointed_properties(::OSIM) = [:clock]
default_checkpointed_properties(::OSIM)  = [:clock]

reset!(model::OSIM) = reset!(model.ocean)

# Make sure to initialize the exchanger here
function initialization_update_state!(model::OSIM)
    initialize!(model.interfaces.exchanger, model.atmosphere)
    update_state!(model)
    return nothing
end

function initialize!(model::OSIM)
    initialize!(model.ocean)
    initialize!(model.interfaces.exchanger, model.atmosphere)
    return nothing
end

checkpointer_address(::OceanSeaIceModel) = "OceanSeaIceModel"

# initialize_jld2_file!(filepath, init, jld2_kw, including, outputs, model::OSIM) =
#     initialize_jld2_file!(filepath, init, jld2_kw, including, outputs, model.ocean.model)

# for prescribed atmosphere just checkpoint the ocean model
# write_output!(c::Checkpointer, model::OSIMPA) = write_output!(c, model.ocean.model)

# function write_output!(c::Checkpointer, model::OSIMPA)
#     atmosphere = model.atmosphere
#     ocean = model.ocean.model

#     write_output!(c, model) # just saves the clock

#     # deals with model components
#     write_output!(c, atmosphere)
#     write_output!(c, ocean)

#     return nothing
# end

function set_clock!(model::OSIM, clock)
    model.clock.time = clock.time
    model.clock.iteration = clock.iteration
    model.clock.last_Δt = clock.last_Δt
    model.clock.last_stage_Δt = clock.last_stage_Δt
    model.clock.stage = clock.stage
    return nothing
end

# function set!(model::OSIM, checkpoint_file_path::AbstractString)
#     atmosphere = model.atmosphere
#     ocean = model.ocean.model

#     set!(model, checkpoint_file_path)

#     # deals with model components
#     set!(ocean, checkpoint_file_path)
#     set!(atmosphere, checkpoint_file_path)

#     return nothing
# end

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

function OceanSeaIceModel(ocean, sea_ice=FreezingLimitedOceanTemperature(eltype(ocean.model));
                          atmosphere = nothing,
                          radiation = Radiation(architecture(ocean.model)),
                          clock = deepcopy(ocean.model.clock),
                          ocean_reference_density = reference_density(ocean),
                          ocean_heat_capacity = heat_capacity(ocean),
                          sea_ice_reference_density = reference_density(sea_ice),
                          sea_ice_heat_capacity = heat_capacity(sea_ice),
                          interfaces = nothing)

    if !isnothing(ocean.callbacks)
        # Remove some potentially irksome callbacks from the ocean simulation
        pop!(ocean.callbacks, :stop_time_exceeded, nothing)
        pop!(ocean.callbacks, :stop_iteration_exceeded, nothing)
        pop!(ocean.callbacks, :wall_time_limit_exceeded, nothing)
        pop!(ocean.callbacks, :nan_checker, nothing)
    end

    if sea_ice isa SeaIceSimulation
        if !isnothing(sea_ice.callbacks)
            pop!(sea_ice.callbacks, :stop_time_exceeded, nothing)
            pop!(sea_ice.callbacks, :stop_iteration_exceeded, nothing)
            pop!(sea_ice.callbacks, :wall_time_limit_exceeded, nothing)
            pop!(sea_ice.callbacks, :nan_checker, nothing)
        end
    end

    # Contains information about flux contributions: bulk formula, prescribed fluxes, etc.
    if isnothing(interfaces) && !(isnothing(atmosphere) && isnothing(sea_ice))
        interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice;
                                         ocean_reference_density,
                                         ocean_heat_capacity,
                                         sea_ice_reference_density,
                                         sea_ice_heat_capacity,
                                         radiation)
    end

    arch = architecture(ocean.model.grid)

    ocean_sea_ice_model = OceanSeaIceModel(arch,
                                           clock,
                                           atmosphere,
                                           sea_ice,
                                           ocean,
                                           interfaces)

    # Make sure the initial temperature of the ocean
    # is not below freezing and above melting near the surface
    above_freezing_ocean_temperature!(ocean, sea_ice)
    initialization_update_state!(ocean_sea_ice_model)

    return ocean_sea_ice_model
end

# Check for NaNs in the first prognostic field (generalizes to prescribed velocities).
function default_nan_checker(model::OSIM)
    u_ocean = model.ocean.model.velocities.u
    nan_checker = NaNChecker((; u_ocean))
    return nan_checker
end

@kernel function _above_freezing_ocean_temperature!(T, grid, S, ℵ, liquidus)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        for k in 1:Nz-1
            Tm = melting_temperature(liquidus, S[i, j, k])
            T[i, j, k] = max(T[i, j, k], Tm)
        end

        ℵi = ℵ[i, j, 1]
        Tm = melting_temperature(liquidus, S[i, j, Nz])
        T[i, j, Nz] = ifelse(ℵi > 0, Tm, T[i, j, Nz])
    end
end

# Fallback
above_freezing_ocean_temperature!(ocean, sea_ice) = nothing

function above_freezing_ocean_temperature!(ocean, sea_ice::SeaIceSimulation)
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S
    ℵ = sea_ice.model.ice_concentration
    liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus

    grid = ocean.model.grid
    arch = architecture(grid)
    launch!(arch, grid, :xy, _above_freezing_ocean_temperature!, T, grid, S, ℵ, liquidus)

    return nothing
end
