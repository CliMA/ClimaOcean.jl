module OceanSeaIceModels

export OceanSeaIceModel, SimilarityTheoryFluxes, FreezingLimitedOceanTemperature
export Radiation, LatitudeDependentAlbedo
export SkinTemperature, BulkTemperature

using Oceananigans
using SeawaterPolynomials

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

function downwelling_radiation end
function freshwater_flux end
function reference_density end
function heat_capacity end

const default_gravitational_acceleration = 9.80665
const default_freshwater_density = 1000

const SeaIceSimulation = Simulation{<:SeaIceModel}

sea_ice_thickness(::Nothing) = ZeroField()
sea_ice_thickness(sea_ice::SeaIceSimulation) = sea_ice.model.ice_thickness

sea_ice_concentration(::Nothing) = ZeroField()
sea_ice_concentration(sea_ice::SeaIceSimulation) = sea_ice.model.ice_concentration

#####
##### Some implementation
#####

# Atmosphere interface
interpolate_atmosphere_state!(interfaces, atmosphere, coupled_model) = nothing
compute_net_atmosphere_fluxes!(coupled_model) = nothing

# TODO: import this last
include("PrescribedAtmospheres.jl")

using .PrescribedAtmospheres:
    PrescribedAtmosphere,
    PrescribedAtmosphereThermodynamicsParameters,
    TwoBandDownwellingRadiation

include("InterfaceComputations/InterfaceComputations.jl")

using .InterfaceComputations

import .InterfaceComputations:
    compute_atmosphere_ocean_fluxes!,
    compute_atmosphere_sea_ice_fluxes!,
    compute_sea_ice_ocean_fluxes!

include("ocean_sea_ice_model.jl")
include("freezing_limited_ocean_temperature.jl")
include("time_step_ocean_sea_ice_model.jl")

# "No atmosphere" implementation
const NoAtmosphereModel = OceanSeaIceModel{<:Any, Nothing}

compute_atmosphere_ocean_fluxes!(::NoAtmosphereModel) = nothing
compute_atmosphere_sea_ice_fluxes!(::NoAtmosphereModel) = nothing

const PrescribedAtmosphereModel = OceanSeaIceModel{<:Any, <:PrescribedAtmosphere}

compute_net_atmosphere_fluxes!(::PrescribedAtmosphereModel) = nothing

# "No sea ice" implementation
const NoSeaIceModel = Union{OceanSeaIceModel{Nothing}, FreezingLimitedCoupledModel}

# Fallback
compute_sea_ice_ocean_fluxes!(::OceanSeaIceModel{Nothing}) = nothing
compute_atmosphere_sea_ice_fluxes!(::NoSeaIceModel) = nothing

# "Only ocean" implementation
const OnlyOceanModel = Union{OceanSeaIceModel{Nothing, Nothing}, OceanSeaIceModel{Nothing, <:FreezingLimitedOceanTemperature}}

compute_atmosphere_sea_ice_fluxes!(::OnlyOceanModel) = nothing
compute_sea_ice_ocean_fluxes!(::OnlyOceanModel) = nothing

end # module
