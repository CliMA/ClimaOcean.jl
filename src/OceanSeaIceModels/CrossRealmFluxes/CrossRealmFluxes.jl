module CrossRealmFluxes

using Oceananigans

export SurfaceRadiation,
       OceanSeaIceSurfaceFluxes

using ..OceanSeaIceModels: SKOFTS

#####
##### Utilities
#####

@inline stateindex(a::Number, i, j, k, time) = a
@inline stateindex(a::SKOFTS, i, j, k, time) = @inbounds a[i, j, k, time]
@inline stateindex(a::AbstractArray, i, j, k, time) = @inbounds a[i, j, k]
@inline Δϕt²(i, j, k, grid, ϕ1, ϕ2, time) = (stateindex(ϕ1, i, j, k, time) - stateindex(ϕ2, i, j, k, time))^2

@inline function stateindex(a::Tuple, i, j, k, time)
    N = length(a)
    ntuple(Val(N)) do n
        stateindex(a[n], i, j, k, time)
    end
end

@inline function stateindex(a::NamedTuple, i, j, k, time)
    vals = stateindex(values(a), i, j, k, time)
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

include("surface_radiation.jl")
include("similarity_theory_surface_fluxes.jl")
include("ocean_sea_ice_surface_fluxes.jl")

# include("ocean_sea_ice_model_fluxes.jl")
# include("ocean_sea_ice_surfaces.jl")
# include("atmosphere_sea_ice_fluxes.jl")
# include("atmosphere_ocean_momentum_flux.jl")
# include("sea_ice_ocean_fluxes.jl")

end # module

