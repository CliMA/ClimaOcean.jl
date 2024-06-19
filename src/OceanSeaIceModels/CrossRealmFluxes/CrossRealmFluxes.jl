module CrossRealmFluxes

using Oceananigans
using Adapt 

export Radiation,
       OceanSeaIceSurfaceFluxes,
       SimilarityTheoryTurbulentFluxes

using ..OceanSeaIceModels: SKOFTS, default_gravitational_acceleration

import ..OceanSeaIceModels: surface_velocities,
                            surface_tracers

#####
##### Utilities
#####

@inline stateindex(a::Number, i, j, k, grid, time) = a
@inline stateindex(a::SKOFTS, i, j, k, grid, time) = @inbounds a[i, j, k, time]
@inline stateindex(a::AbstractArray, i, j, k, grid, time) = @inbounds a[i, j, k]
@inline Δϕt²(i, j, k, grid, ϕ1, ϕ2, time) = (stateindex(ϕ1, i, j, k, grid, time) - stateindex(ϕ2, i, j, k, grid, time))^2

@inline function stateindex(a::Tuple, i, j, k, grid, time)
    N = length(a)
    ntuple(Val(N)) do n
        stateindex(a[n], i, j, k, grid, time)
    end
end

@inline function stateindex(a::NamedTuple, i, j, k, grid, time)
    vals = stateindex(values(a), i, j, k, grid, time)
    names = keys(a)
    return NamedTuple{names}(vals)
end

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
include("sea_ice_ocean_fluxes.jl")

end # module
