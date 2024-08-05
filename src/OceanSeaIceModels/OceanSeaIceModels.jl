module OceanSeaIceModels

export OceanSeaIceModel, SimilarityTheoryTurbulentFluxes, Radiation

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

using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

using ClimaOcean: stateindex

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

function surface_velocities end
function surface_tracers end
function surface_horizontal_velocities end
function surface_active_tracers end
function downwelling_radiation end
function freshwater_flux end
function reference_density end
function heat_capacity end

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

include("minimum_temperature_sea_ice.jl")
include("ocean_sea_ice_model.jl")
include("ocean_only_model.jl")
include("time_step_ocean_sea_ice_model.jl")

import .CrossRealmFluxes:
    compute_atmosphere_ocean_fluxes!,
    compute_sea_ice_ocean_fluxes!

# "No atmosphere" implementation
const NoAtmosphereModel = OceanSeaIceModel{<:Any, Nothing}

compute_atmosphere_ocean_fluxes!(coupled_model::NoAtmosphereModel) = nothing

sea_ice_thickness(sea_ice::SeaIceModel)     = sea_ice.model.ice_thickness.data
sea_ice_concentration(sea_ice::SeaIceModel) = sea_ice.ice_concentration.data
sea_ice_salinity(sea_ice::SeaIceModel)      = sea_ice.tracers.S.data
sea_ice_velocities(sea_ice::SeaIceModel)    = map(u -> u.data, sea_ice.velocities)

sea_ice_thickness(::Nothing) = nothing
sea_ice_concentration(::Nothing) = nothing
sea_ice_salinity(::Nothing) = nothing
sea_ice_velocities(::Nothing) = nothing
sea_ice_temperature()
end # module

