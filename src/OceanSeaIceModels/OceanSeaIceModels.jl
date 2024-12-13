module OceanSeaIceModels

export OceanSeaIceModel, SimilarityTheoryFluxes, FreezingLimitedOceanTemperature
export Radiation, LatitudeDependentAlbedo
export SkinTemperature, BulkTemperature

using Oceananigans
using SeawaterPolynomials

using Oceananigans.Operators

using Oceananigans.Utils: launch!, Time
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!, BoundaryCondition
using Oceananigans.Grids: architecture
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

sea_ice_thickness(::Nothing) = nothing
sea_ice_thickness(sea_ice::SeaIceSimulation) = sea_ice.model.ice_thickness

sea_ice_concentration(::Nothing) = nothing
sea_ice_concentration(sea_ice::SeaIceSimulation) = sea_ice.model.ice_concentration

#####
##### Some implementation
#####

include("PrescribedAtmospheres.jl")

using .PrescribedAtmospheres:
    PrescribedAtmosphere,
    PrescribedAtmosphereThermodynamicsParameters,
    TwoBandDownwellingRadiation

include("CrossRealmFluxes/CrossRealmFluxes.jl")

using .CrossRealmFluxes

import .CrossRealmFluxes:
    compute_atmosphere_ocean_fluxes!,
    compute_atmosphere_sea_ice_fluxes!,
    compute_sea_ice_ocean_fluxes!,
    limit_fluxes_over_sea_ice!

include("ocean_sea_ice_model.jl")
include("freezing_limited_ocean_temperature.jl")
include("time_step_ocean_sea_ice_model.jl")

# "No atmosphere" implementation
const NoAtmosphereModel = OceanSeaIceModel{<:Any, Nothing}

compute_atmosphere_ocean_fluxes!(::NoAtmosphereModel) = nothing
compute_atmosphere_sea_ice_fluxes!(::NoAtmosphereModel) = nothing

const NoSeaIceModel = {OceanSeaIceModel{Nothing}, FreezingLimitedCoupledModel}

# Fallback
compute_sea_ice_ocean_fluxes!(::NoSeaIceModel) = nothing
compute_atmosphere_ocean_fluxes!(::NoSeaIceModel) = nothing

end # module
