module IceOceanModel

using Oceananigans.Operators

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
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

const ℒₑ = 2.5e6 # J/kg Latent heat of evaporation
const σᴮ = 5.67e-8 # W/m²/K⁴ Stefan-Boltzmann constant

include("ice_ocean_model.jl")
include("ice_ocean_atmosphere_fluxes.jl")
include("only_ocean_model_fluxes.jl")
include("AtmosphericForcings.jl")

using .AtmosphericForcings

# Check for NaNs in the first prognostic field (generalizes to prescribed velocitries).
function default_nan_checker(model::IceOceanModel)
    u_ocean = model.ocean.model.velocities.u
    nan_checker = NaNChecker((; u_ocean))
    return nan_checker
end

end # module
