module OceanSeaIceModels

using Oceananigans.Operators

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!, BoundaryCondition
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: tick!
using Oceananigans.Utils: launch!

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

# We should not declare these; they need to be settable.
# const ℒₑ = 2.5e6 # J/kg Latent heat of evaporation
# const σᴮ = 5.67e-8 # W/m²/K⁴ Stefan-Boltzmann constant

using Oceananigans
using Oceananigans.Utils: Time
using Oceananigans.Grids: architecture
using Oceananigans.Models: AbstractModel

#####
##### Interface
#####

function surface_velocities(ocean::Simulation{<:HydrostaticFreeSurfaceModel})
    grid = ocean.model.grid
    Nz = size(grid, 3)
    u = interior(ocean.model.velocities.u, :, :, Nz)
    v = interior(ocean.model.velocities.v, :, :, Nz)
    w = interior(ocean.model.velocities.w, :, :, Nz+1)
    return (; u, v, w)
end

function surface_flux(f::Field)
    top_bc = f.boundary_conditions.top
    if top_bc isa BoundaryCondition{<:Oceananigans.BoundaryConditions.Flux}
        return top_bc.condition
    else
        return nothing
    end
end

#####
##### Some implementation
#####

include("atmosphere_sea_ice_fluxes.jl")
include("atmosphere_ocean_fluxes.jl")
include("sea_ice_ocean_fluxes.jl")
include("ocean_sea_ice_model.jl")
include("ocean_only_model.jl")

# "No atmosphere" implementation
const NoAtmosphereModel = OceanSeaIceModel{<:Any, <:Any, Nothing}
compute_atmosphere_ocean_fluxes!(coupled_model::NoAtmosphereModel) = nothing

include("PrescribedAtmospheres.jl")

using .PrescribedAtmospheres: PrescribedAtmosphere

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

