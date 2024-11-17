module OceanSeaIceModels

export OceanSeaIceModel, SimilarityTheoryTurbulentFluxes, FreezingLimitedOceanTemperature
export Radiation, LatitudeDependentAlbedo

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

using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

using ClimaOcean: stateindex

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

function downwelling_radiation end
function freshwater_flux end
function reference_density end
function heat_capacity end

sea_ice_thickness(::Nothing) = nothing
sea_ice_concentration(::Nothing) = nothing

const default_gravitational_acceleration = 9.80665
const default_freshwater_density = 1000

#####
##### Some implementation
#####

include("PrescribedAtmospheres.jl")

using .PrescribedAtmospheres:
    PrescribedAtmosphere,
    TwoBandDownwellingRadiation

include("CrossRealmFluxes/CrossRealmFluxes.jl")

using .CrossRealmFluxes

import .CrossRealmFluxes:
    compute_atmosphere_ocean_fluxes!,
    compute_sea_ice_ocean_fluxes!,
    limit_fluxes_over_sea_ice!

include("ocean_sea_ice_model.jl")
include("freezing_limited_ocean_temperature.jl")
include("time_step_ocean_sea_ice_model.jl")

# "No atmosphere" implementation
const NoAtmosphereModel = OceanSeaIceModel{<:Any, Nothing}

compute_atmosphere_ocean_fluxes!(coupled_model::NoAtmosphereModel) = nothing

end # module
