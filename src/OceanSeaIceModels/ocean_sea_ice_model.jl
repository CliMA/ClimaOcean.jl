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
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!, time
import Oceananigans.Utils: prettytime
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker, initialization_update_state!

mutable struct OceanSeaIceModel{I, A, O, F, C, Arch} <: AbstractModel{Nothing, Arch}
    architecture :: Arch
    clock :: C
    atmosphere :: A
    sea_ice :: I
    ocean :: O
    interfaces :: F
end

const OSIM = OceanSeaIceModel

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
architecture(model::OSIM)           = architecture(model.ocean.model)
Base.eltype(model::OSIM)            = Base.eltype(model.ocean.model)
prettytime(model::OSIM)             = prettytime(model.clock.time)
iteration(model::OSIM)              = model.clock.iteration
timestepper(::OSIM)                 = nothing
default_included_properties(::OSIM) = tuple()
prognostic_fields(cm::OSIM)         = nothing
fields(::OSIM)                      = NamedTuple()
default_clock(TT)                   = Oceananigans.TimeSteppers.Clock{TT}(0, 0, 1)

function reset!(model::OSIM)
    reset!(model.ocean)
    return nothing
end

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

"""
    OceanSeaIceModel(ocean, sea_ice=FreezingLimitedOceanTemperature(eltype(ocean.model));
                     atmosphere = nothing,
                     radiation = Radiation(architecture(ocean.model)),
                     clock = deepcopy(ocean.model.clock),
                     ocean_reference_density = reference_density(ocean),
                     ocean_heat_capacity = heat_capacity(ocean),
                     sea_ice_reference_density = reference_density(sea_ice),
                     sea_ice_heat_capacity = heat_capacity(sea_ice),
                     interfaces = nothing)

Construct a coupled ocean-sea ice model that simulates the interaction between ocean and sea ice components.

# Arguments
- `ocean`: A representation of a possibly time-dependent ocean state. Currently, only `Oceananigans.Simulation`s
           of `Oceananigans.HydrostaticFreeSurfaceModel` are tested.
- `sea_ice`: A representation of a possibly time-dependent sea ice state.
             For example, the minimalist `FreezingLimitedOceanTemperature` represents
             oceanic latent heating during freezing only, but does not evolve sea ice variables.
             For prognostica sea ice use an `Oceananigans.Simulation`s of `ClimaSeaIce.SeaIceModel`.

# Keyword Arguments
- `atmosphere`: A representation of a possibly time-dependent atmospheric state. Default: `nothing`.
- `radiation`: Radiation component used to compute surface fluxes at the bottom of the atmosphere.
- `clock`: Keeps track of time.
- `ocean_reference_density`: Reference density for the ocean. Defaults to value from ocean model
- `ocean_heat_capacity`: Heat capacity for the ocean. Defaults to value from ocean model
- `sea_ice_reference_density`: Reference density for sea ice. Defaults to value from sea ice model
- `sea_ice_heat_capacity`: Heat capacity for sea ice. Defaults to value from sea ice model
- `interfaces`: Component interfaces for coupling. Defaults to `nothing` and will be constructed automatically

# Stability Functions
The model uses similarity theory for turbulent fluxes between components. You can customize the stability functions
by creating a new `SimilarityTheoryFluxes` object with your desired stability functions. For example:

```jldoctest ocean_sea_ice_model
using ClimaOcean
using Oceananigans

grid = RectilinearGrid(size=10, z=(-100, 0), topology=(Flat, Flat, Bounded))
ocean = ocean_simulation(grid)

# Three choices for stability function:
# "No stability function", which also apply to neutral boundary layers
stability_functions = nothing

# Atmosphere-sea ice specific stability functions
stability_functions = ClimaOcean.OceanSeaIceModels.atmosphere_sea_ice_stability_functions(Float64)

# Edson et al. (2013) stability functions
stability_functions = ClimaOcean.OceanSeaIceModels.atmosphere_ocean_stability_functions(Float64)

atmosphere_ocean_flux_formulation = SimilarityTheoryFluxes(; stability_functions)
interfaces = ClimaOcean.OceanSeaIceModels.ComponentInterfaces(nothing, ocean; atmosphere_ocean_flux_formulation)
model = OceanSeaIceModel(ocean; interfaces)

# output
OceanSeaIceModel{CPU}(time = 0 seconds, iteration = 0)
├── ocean: HydrostaticFreeSurfaceModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── atmosphere: Nothing
├── sea_ice: FreezingLimitedOceanTemperature{ClimaSeaIce.SeaIceThermodynamics.LinearLiquidus{Float64}}
└── interface: ComponentInterfaces
```

The available stability function options include:
- `atmosphere_ocean_stability_functions`: Based on Edson et al. (2013)
- `atmosphere_sea_ice_stability_functions`: Specifically designed for atmosphere-sea ice interactions
- `nothing`: No stability functions will be used
- Custom stability functions can be created by defining functions of the "stability parameter" 
  (the flux Richardson number), `ζ`.
```
"""
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

time(coupled_model::OceanSeaIceModel) = coupled_model.clock.time

# Check for NaNs in the first prognostic field (generalizes to prescribed velocities).
function default_nan_checker(model::OceanSeaIceModel)
    u_ocean = model.ocean.model.velocities.u
    nan_checker = NaNChecker((; u_ocean))
    return nan_checker
end

@kernel function _above_freezing_ocean_temperature!(T, grid, S, ℵ, liquidus)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        for k in 1:Nz
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
