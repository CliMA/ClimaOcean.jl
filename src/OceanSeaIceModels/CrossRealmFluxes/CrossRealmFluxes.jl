module CrossRealmFluxes

using Oceananigans
using Adapt 

export Radiation,
       OceanSeaIceSurfaceFluxes,
       SimilarityTheoryTurbulentFluxes

using ..OceanSeaIceModels: default_gravitational_acceleration

import ..OceanSeaIceModels: surface_velocities,
                            surface_tracers

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

function surface_velocities(ocean::Simulation{<:HydrostaticFreeSurfaceModel})
    grid = ocean.model.grid
    Nz = size(grid, 3)
    u = view(ocean.model.velocities.u.data, :, :, Nz)
    v = view(ocean.model.velocities.v.data, :, :, Nz)
    w = view(ocean.model.velocities.w.data, :, :, Nz+1)
    return (; u, v, w)
end

function surface_tracers(ocean::Simulation{<:HydrostaticFreeSurfaceModel})
    grid = ocean.model.grid
    Nz = size(grid, 3)
    tracers = ocean.model.tracers
    names = keys(tracers)
    sfc_tracers = NamedTuple(name => view(tracers[name].data, :, :, Nz) for name in names)
    return sfc_tracers
end

include("radiation.jl")
include("tabulated_albedo.jl")
include("roughness_lengths.jl")
include("stability_functions.jl")
include("seawater_saturation_specific_humidity.jl")
include("similarity_theory_turbulent_fluxes.jl")
include("ocean_sea_ice_surface_fluxes.jl")
include("atmosphere_ocean_fluxes.jl")
include("atmosphere_sea_ice_fluxes.jl")
include("sea_ice_ocean_fluxes.jl")

end # module
