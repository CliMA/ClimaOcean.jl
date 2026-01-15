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
    compute_sea_ice_ocean_fluxes!,
    # Sea ice-ocean heat flux formulations
    IceBathHeatFlux,
    ThreeEquationHeatFlux,
    # Friction velocity formulations
    MomentumBasedFrictionVelocity

using ..OceanSeaIceModels: default_gravitational_acceleration,
                           default_freshwater_density,
                           thermodynamics_parameters,
                           surface_layer_height,
                           boundary_layer_height

import ClimaOcean: stateindex
import Oceananigans.Simulations: initialize!

#####
##### Functions extended by component models
#####

net_fluxes(::Nothing) = nothing

#####
##### Utilities
#####

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

# Sea ice-ocean heat flux formulations
include("friction_velocity.jl")
include("sea_ice_ocean_heat_flux_formulations.jl")

include("component_interfaces.jl")
include("atmosphere_ocean_fluxes.jl")
include("atmosphere_sea_ice_fluxes.jl")
include("sea_ice_ocean_fluxes.jl")

end # module
