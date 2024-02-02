module OceanSeaIceModels

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

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

const SomeKindOfFieldTimeSeries = Union{FieldTimeSeries,
                                        GPUAdaptedFieldTimeSeries}

const SKOFTS = SomeKindOfFieldTimeSeries

function surface_velocities end
function surface_tracers end
function downwelling_radiation end
function freshwater_flux end
function reference_density end
function heat_capacity end

sea_ice_thickness(::Nothing) = nothing
sea_ice_concentration(::Nothing) = nothing

#####
##### Some implementation
#####

include("PrescribedAtmospheres.jl")

using .PrescribedAtmospheres:
    PrescribedAtmosphere,
    TwoStreamDownwellingRadiation

include("CrossRealmFluxes/CrossRealmFluxes.jl")

using .CrossRealmFluxes

include("ocean_sea_ice_model.jl")
include("ocean_only_model.jl")
include("time_step_ocean_sea_ice_model.jl")

import .CrossRealmFluxes: compute_atmosphere_ocean_fluxes!

# "No atmosphere" implementation
const NoAtmosphereModel = OceanSeaIceModel{<:Any, Nothing}
compute_atmosphere_ocean_fluxes!(coupled_model::NoAtmosphereModel) = nothing

end # module

