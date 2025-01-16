module CrossRealmFluxes

using Oceananigans
using Oceananigans.Fields: AbstractField
using Adapt 

export Radiation,
       CrossRealmSurfaceFluxes,
       LatitudeDependentAlbedo,
       SimilarityTheoryFluxes,
       SkinTemperature, 
       BulkTemperature

using ..OceanSeaIceModels: default_gravitational_acceleration

import ClimaOcean: stateindex

#####
##### Utilities
#####

const c = Center()
const f = Face()

function surface_computations_kernel_parameters(grid)
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
include("stability_functions.jl")
# include("seawater_saturation_specific_humidity.jl")
include("interface_states.jl")
include("similarity_theory_turbulent_fluxes.jl")
include("constant_coefficient_turbulent_fluxes.jl")

# Total fluxes
include("cross_realm_surface_fluxes.jl")
include("interpolate_atmospheric_state.jl")
include("assemble_interface_fluxes.jl")
include("atmosphere_ocean_fluxes.jl")
include("atmosphere_sea_ice_fluxes.jl")
include("sea_ice_ocean_fluxes.jl")

end # module
