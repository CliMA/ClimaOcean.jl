module CrossRealmFluxes

using Oceananigans
using Adapt 

export Radiation,
       OceanSeaIceSurfaceFluxes,
       LatitudeDependentAlbedo,
       SimilarityTheoryTurbulentFluxes,
       SkinTemperature, 
       BulkTemperature

using ..OceanSeaIceModels: default_gravitational_acceleration

import ClimaOcean: stateindex

#####
##### Utilities
#####

function surface_flux(f::Field)
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
include("state_differences.jl")
include("roughness_lengths.jl")
include("stability_functions.jl")
include("seawater_saturation_specific_humidity.jl")
include("surface_temperature.jl")
include("similarity_theory_turbulent_fluxes.jl")
include("constant_coefficient_turbulent_fluxes.jl")

# Total fluxes
include("ocean_sea_ice_surface_fluxes.jl")
include("compute_atmosphere_surface_fluxes.jl")
include("atmosphere_ocean_fluxes.jl")
include("atmosphere_sea_ice_fluxes.jl")
include("sea_ice_ocean_fluxes.jl")

end # module
