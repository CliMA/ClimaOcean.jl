module InterfaceComputations

using Oceananigans
using Oceananigans.Fields: AbstractField
using Oceananigans.Utils: KernelParameters
using Adapt

export
    Radiation,
    ComponentInterfaces,
    LatitudeDependentAlbedo,
    SimilarityTheoryFluxes,
    MomentumRoughnessLength,
    ScalarRoughnessLength,
    CoefficientBasedFluxes,
    SkinTemperature,
    BulkTemperature,
    atmosphere_ocean_stability_functions,
    atmosphere_sea_ice_stability_functions,
    compute_atmosphere_ocean_fluxes!,
    compute_atmosphere_sea_ice_fluxes!,
    compute_sea_ice_ocean_flu

using ..OceanSeaIceModels: default_gravitational_acceleration,
                           default_freshwater_density

import ClimaOcean: stateindex

import Oceananigans.Simulations: initialize!

import ..OceanSeaIceModels:
    compute_net_atmosphere_fluxes!,
    compute_net_sea_ice_fluxes!,
    compute_net_ocean_fluxes!

#####
##### Utilities
#####

const c = Center()
const f = Face()

function interface_kernel_parameters(grid)
    Nx, Ny, Nz = size(grid)
    single_column_grid = Nx == 1 && Ny == 1

    if single_column_grid
        kernel_parameters = KernelParameters(1:1, 1:1)
    else
        # Compute fluxes into halo regions, ie from 0:Nx+1 and 0:Ny+1
        kernel_parameters = KernelParameters(0:Nx+1, 0:Ny+1)
    end

    return kernel_parameters
end

# Needs to be extended by each component model
net_fluxes(::Nothing) = nothing

function surface_flux(f::AbstractField)
    top_bc = f.boundary_conditions.top
    if top_bc isa BoundaryCondition{<:Oceananigans.BoundaryConditions.Flux}
        return top_bc.condition
    else
        return nothing
    end
end

# Radiation
include("radiation.jl")
include("latitude_dependent_albedo.jl")
include("tabulated_albedo.jl")

# Turbulent fluxes
include("roughness_lengths.jl")
include("interface_states.jl")
include("compute_interface_state.jl")
include("similarity_theory_turbulent_fluxes.jl")
include("coefficient_based_turbulent_fluxes.jl")

# State exchanger and interfaces
include("state_exchanger.jl")
include("component_interfaces.jl")
include("atmosphere_ocean_fluxes.jl")
include("atmosphere_sea_ice_fluxes.jl")
include("sea_ice_ocean_fluxes.jl")

end # module
