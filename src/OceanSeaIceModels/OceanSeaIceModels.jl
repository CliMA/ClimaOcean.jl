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
    compute_sea_ice_ocean_fluxes!

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
import Oceananigans: fields, prognostic_fields
import Oceananigans.Architectures: architecture
import Oceananigans.Fields: set!
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker, initialization_update_state!
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
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
#####  Fallbacks for a single-component model
#####

#                                      Sea Ice   | Atmosphere | Ocean
const OnlyOceanModel      = OceanSeaIceModel{<:Nothing, <:Nothing, <:Any}
const OnlyAtmosphereModel = OceanSeaIceModel{<:Nothing, <:Any,     <:Nothing}
const OnlySeaIceModel     = OceanSeaIceModel{<:Any,     <:Nothing, <:Nothing}

SingleComponentModel = Union{OnlyOceanModel, OnlyAtmosphereModel, OnlySeaIceModel}

InterfaceComputations.compute_atmosphere_ocean_fluxes!(::SingleComponentModel) = nothing
InterfaceComputations.compute_atmosphere_sea_ice_fluxes!(::SingleComponentModel) = nothing
InterfaceComputations.compute_sea_ice_ocean_fluxes!(::SingleComponentModel) = nothing

end # module
