module OceanSeaIceModels

export
    OceanSeaIceModel,
    SimilarityTheoryFluxes,
    CoefficientBasedFluxes,
    FreezingLimitedOceanTemperature,
    Radiation,
    LatitudeDependentAlbedo,
    SkinTemperature,
    BulkTemperature,
    compute_atmosphere_ocean_fluxes!,
    compute_atmosphere_sea_ice_fluxes!,
    compute_sea_ice_ocean_fluxes!,
    # Sea ice-ocean heat flux formulations
    IceBathHeatFlux,
    ThreeEquationHeatFlux,
    # Friction velocity formulations
    MomentumBasedFrictionVelocity

using Oceananigans
using Oceananigans.Operators
using Oceananigans.Utils: launch!, KernelParameters
using Oceananigans.Units: Time
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!, BoundaryCondition
using Oceananigans.Grids: architecture
using Oceananigans.Fields: ZeroField
using Oceananigans.TimeSteppers: tick!
using Oceananigans.Models: AbstractModel
using Oceananigans.OutputReaders: FieldTimeSeries, GPUAdaptedFieldTimeSeries

using ClimaSeaIce: SeaIceModel
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

using ClimaOcean: stateindex

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

import Thermodynamics as AtmosphericThermodynamics

# Simulations interface
import Oceananigans: fields, prognostic_fields, prognostic_state, restore_prognostic_state!
import Oceananigans.Architectures: architecture
import Oceananigans.Fields: set!
import Oceananigans.Models: NaNChecker, default_nan_checker, initialization_update_state!
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: timestepper, reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!, time
import Oceananigans.Utils: prettytime

include("components.jl")

#####
##### The coupled model
#####

const default_gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration
const default_freshwater_density = 1000 # kg m⁻³

include("InterfaceComputations/InterfaceComputations.jl")

using .InterfaceComputations

include("ocean_sea_ice_model.jl")
include("time_step_ocean_sea_ice_model.jl")

#####
#####  Fallbacks for no-interface models
#####

using .InterfaceComputations: ComponentInterfaces, AtmosphereInterface, SeaIceOceanInterface

const NoSeaIceInterface = ComponentInterfaces{<:AtmosphereInterface,  <:Nothing, <:Nothing}
const NoOceanInterface  = ComponentInterfaces{<:Nothing, <:AtmosphereInterface,  <:Nothing}
const NoAtmosInterface  = ComponentInterfaces{<:Nothing, <:Nothing, <:SeaIceOceanInterface}
const NoInterface       = ComponentInterfaces{<:Nothing, <:Nothing, <:Nothing}

const NoSeaIceInterfaceModel = OceanSeaIceModel{I, A, O, <:NoSeaIceInterface} where {I, A, O}
const NoAtmosInterfaceModel  = OceanSeaIceModel{I, A, O, <:NoAtmosInterface}  where {I, A, O}
const NoOceanInterfaceModel  = OceanSeaIceModel{I, A, O, <:NoOceanInterface}  where {I, A, O}
const NoInterfaceModel       = OceanSeaIceModel{I, A, O, <:NoInterface}  where {I, A, O}

InterfaceComputations.compute_atmosphere_sea_ice_fluxes!(::NoSeaIceInterfaceModel) = nothing
InterfaceComputations.compute_sea_ice_ocean_fluxes!(::NoSeaIceInterfaceModel) = nothing

InterfaceComputations.compute_atmosphere_ocean_fluxes!(::NoAtmosInterfaceModel) = nothing
InterfaceComputations.compute_atmosphere_sea_ice_fluxes!(::NoAtmosInterfaceModel) = nothing

InterfaceComputations.compute_atmosphere_ocean_fluxes!(::NoOceanInterfaceModel) = nothing
InterfaceComputations.compute_sea_ice_ocean_fluxes!(::NoOceanInterfaceModel) = nothing

InterfaceComputations.compute_atmosphere_ocean_fluxes!(::NoInterfaceModel) = nothing
InterfaceComputations.compute_atmosphere_sea_ice_fluxes!(::NoInterfaceModel) = nothing
InterfaceComputations.compute_sea_ice_ocean_fluxes!(::NoInterfaceModel) = nothing

end # module
