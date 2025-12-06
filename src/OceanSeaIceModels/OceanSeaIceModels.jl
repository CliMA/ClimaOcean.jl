module OceanSeaIceModels

export
    OceanSeaIceModel,
    SimilarityTheoryFluxes,
    CoefficientBasedFluxes,
    FreezingLimitedOceanTemperature,
    Radiation,
    LatitudeDependentAlbedo,
    SkinTemperature,
    BulkTemperature

using SeawaterPolynomials

using Oceananigans
using Oceananigans.Operators
using Oceananigans.Utils: launch!, Time, KernelParameters
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

function downwelling_radiation end
function freshwater_flux end
function reference_density end
function heat_capacity end

const default_gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration
const default_freshwater_density = 1000 # kg m⁻³

# Our default ocean and sea ice models
const SeaIceSimulation = Simulation{<:SeaIceModel}
const OceananigansSimulation = Simulation{<:HydrostaticFreeSurfaceModel}

Base.eltype(ocean::OceananigansSimulation) = eltype(ocean.model.grid)

sea_ice_thickness(::Nothing) = ZeroField()
sea_ice_thickness(sea_ice::SeaIceSimulation) = sea_ice.model.ice_thickness

sea_ice_concentration(::Nothing) = ZeroField()
sea_ice_concentration(sea_ice::SeaIceSimulation) = sea_ice.model.ice_concentration

mutable struct OceanSeaIceModel{I, A, O, F, C, Arch} <: AbstractModel{Nothing, Arch}
    architecture :: Arch
    clock :: C
    atmosphere :: A
    sea_ice :: I
    ocean :: O
    interfaces :: F
end

struct FreezingLimitedOceanTemperature{L}
    liquidus :: L
end

const OSIM = OceanSeaIceModel
const NoAtmosphereModel = OceanSeaIceModel{<:Any, Nothing}
const NoSeaIceModel = Union{OceanSeaIceModel{Nothing}, OceanSeaIceModel{<:FreezingLimitedOceanTemperature}}

#####
##### Some implementation
#####

# Atmosphere interface
interpolate_atmosphere_state!(interfaces, atmosphere, coupled_model) = nothing

# Compute net fluxes:
compute_net_sea_ice_fluxes!(coupled_model,    ::Nothing) = nothing
compute_net_ocean_fluxes!(coupled_model,      ::Nothing) = nothing
compute_net_atmosphere_fluxes!(coupled_model, ::Nothing) = nothing

# "No atmosphere" implementation
compute_atmosphere_ocean_fluxes!(::NoAtmosphereModel) = nothing
compute_atmosphere_sea_ice_fluxes!(::NoAtmosphereModel) = nothing

# "No sea ice" implementation
compute_sea_ice_ocean_fluxes!(::OceanSeaIceModel{Nothing}) = nothing
compute_atmosphere_sea_ice_fluxes!(::NoSeaIceModel) = nothing

# "Only ocean" implementation
const OnlyOceanModel = Union{OceanSeaIceModel{Nothing, Nothing}, OceanSeaIceModel{<:FreezingLimitedOceanTemperature, Nothing}}

compute_atmosphere_sea_ice_fluxes!(::OnlyOceanModel) = nothing
compute_sea_ice_ocean_fluxes!(::OnlyOceanModel) = nothing
compute_net_ocean_fluxes!(::OnlyOceanModel, ocean) = nothing

include("freezing_limited_ocean_temperature.jl")

# TODO: import this last
include("PrescribedAtmospheres.jl")

using .PrescribedAtmospheres:
    PrescribedAtmosphere,
    AtmosphereThermodynamicsParameters,
    TwoBandDownwellingRadiation

include("InterfaceComputations/InterfaceComputations.jl")

using .InterfaceComputations

include("ocean_sea_ice_model.jl")
include("time_step_ocean_sea_ice_model.jl")

end # module
