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

#####
##### Functions extended by sea-ice and ocean models
#####

reference_density(::Nothing) = 0
heat_capacity(::Nothing) = 0

#####
##### Functions extended by sea-ice models
#####

sea_ice_thickness(::Nothing) = ZeroField()
sea_ice_concentration(::Nothing) = ZeroField()

#####
##### Functions extended by atmosphere models
#####

function downwelling_radiation end
function freshwater_flux end
function thermodynamics_parameters end
function surface_layer_height end
function boundary_layer_height end

#####
##### Functions extended by all component models
#####

function interpolate_state! end
function compute_net_fluxes! end

#####
##### The coupled model
##### 

const default_gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration
const default_freshwater_density = 1000 # kg m⁻³

mutable struct OceanSeaIceModel{I, A, O, F, C, Arch} <: AbstractModel{Nothing, Arch}
    architecture :: Arch
    clock :: C
    atmosphere :: A
    sea_ice :: I
    ocean :: O
    interfaces :: F
end

const OSIM = OceanSeaIceModel
const NoAtmosphereModel = OceanSeaIceModel{<:Any, Nothing}

#####
##### Some implementation
#####

# Compute net fluxes:
compute_net_sea_ice_fluxes!(coupled_model,    ::Nothing) = nothing
compute_net_ocean_fluxes!(coupled_model,      ::Nothing) = nothing
compute_net_atmosphere_fluxes!(coupled_model, ::Nothing) = nothing

include("InterfaceComputations/InterfaceComputations.jl")

using .InterfaceComputations

include("ocean_sea_ice_model.jl")
include("time_step_ocean_sea_ice_model.jl")

end # module
