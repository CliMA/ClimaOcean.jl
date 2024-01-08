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

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.Models: timestepper, NaNChecker, default_nan_checker
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!, time
import Oceananigans.Utils: prettytime

const SomeKindOfFieldTimeSeries = Union{FieldTimeSeries,
                                        GPUAdaptedFieldTimeSeries}

const SKOFTS = SomeKindOfFieldTimeSeries

function surface_velocities end
function surface_tracers end
function sea_ice_thickness end
function downwelling_radiation end
function freshwater_flux end
function heat_capacity end
function density end

#####
##### Some implementation
#####

include("CrossRealmFluxes/CrossRealmFluxes.jl")

using .CrossRealmFluxes

include("compute_atmosphere_ocean_fluxes.jl")
include("ocean_sea_ice_model.jl")
include("ocean_only_model.jl")

# "No atmosphere" implementation
const NoAtmosphereModel = OceanSeaIceModel{<:Any, <:Any, Nothing}
compute_atmosphere_ocean_fluxes!(coupled_model::NoAtmosphereModel) = nothing

include("PrescribedAtmospheres.jl")

using .PrescribedAtmospheres: PrescribedAtmosphere, TwoStreamDownwellingRadiation

# Or "AtmosphereModels"
# include("Atmospheres.jl")
# using .Atmospheres

# Check for NaNs in the first prognostic field (generalizes to prescribed velocitries).
function default_nan_checker(model::OceanSeaIceModel)
    u_ocean = model.ocean.model.velocities.u
    nan_checker = NaNChecker((; u_ocean))
    return nan_checker
end

end # module

