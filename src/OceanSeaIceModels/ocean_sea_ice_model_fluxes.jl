using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel
using ClimaSeaIce.SlabSeaIceModels: SlabSeaIceModel

#####
##### Utilities
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
##### Convenience containers for surface fluxes
##### 
##### "Cross realm fluxes" can refer to the flux _data_ (ie, fields representing
##### the total flux for a given variable), or to the flux _components_ / formula.
#####

struct CrossRealmFluxes{M, H, T}
    momentum :: M
    heat :: H
    tracers :: T
end

CrossRealmFluxes(; momentum=nothing, heat=nothing, tracers=nothing) =
    CrossRealmFluxes(momentum, heat, tracers)

Base.summary(osf::CrossRealmFluxes) = "CrossRealmFluxes"
Base.show(io::IO, osf::CrossRealmFluxes) = print(io, summary(osf))

#####
##### Extractors for differnet models (maybe this belongs in the model repo's)
#####

function extract_top_surface_fluxes(model::HydrostaticFreeSurfaceModel)
    u_flux = surface_flux(model.velocities.u)
    v_flux = surface_flux(model.velocities.v)

    ocean_momentum_fluxes = (u = u_flux, v = v_flux)

    ocean_tracers = model.tracers

    ocean_tracer_fluxes = NamedTuple(name => surface_flux(ocean_tracers[name])
                                     for name in keys(ocean_tracers)
                                     if surface_flux(ocean_tracers[name]) isa AbstractArray)

    ocean_fluxes = CrossRealmFluxes(momentum = ocean_momentum_fluxes,
                                    tracers = ocean_tracer_fluxes)

    return ocean_fluxes
end

extract_top_surface_fluxes(model::SlabSeaIceModel) = nothing
extract_bottom_surface_fluxes(model::SlabSeaIceModel) = nothing

#####
##### Total flux across each surface
#####

struct OceanSeaIceSurfaces{O, IT, IB}
    ocean :: O
    sea_ice_top :: IT
    sea_ice_bottom :: IB
end

Base.summary(osis::OceanSeaIceSurfaces) = "OceanSeaIceSurfaces"
Base.show(io::IO, osis::OceanSeaIceSurfaces) = print(io, summary(osis))

function OceanSeaIceSurfaces(ocean, sea_ice=nothing)
    ocean_fluxes = extract_top_surface_fluxes(ocean.model)
            
    if isnothing(sea_ice)
        sea_ice_top_fluxes = nothing
        sea_ice_bottom_fluxes = nothing
    else
        sea_ice_top_fluxes = extract_top_surface_fluxes(sea_ice.model)
        sea_ice_bottom_fluxes = extract_bottom_surface_fluxes(sea_ice.model)
    end

    return OceanSeaIceSurfaces(ocean_fluxes,
                               sea_ice_top_fluxes,  
                               sea_ice_bottom_fluxes)
end

#####
##### Container for organizing information related to fluxes
#####

struct OceanSeaIceModelFluxes{S, R, AO, AI, IO}
    surfaces :: S
    radiation :: R
    atmosphere_ocean :: AO
    atmosphere_sea_ice :: AI
    sea_ice_ocean :: IO
end

function OceanSeaIceModelFluxes(ocean, sea_ice=nothing;
                                radiation = nothing,
                                atmosphere_ocean = nothing,
                                atmosphere_sea_ice = nothing,
                                sea_ice_ocean = nothing)

    surfaces = OceanSeaIceSurfaces(ocean, sea_ice)

    if isnothing(atmosphere_ocean) # defaults
        FT = eltype(ocean.model.grid)
        τˣ = BulkFormula(FT, transfer_coefficient=1e-3)
        τʸ = BulkFormula(FT, transfer_coefficient=1e-3)
        momentum_flux_formulae = (u=τˣ, v=τʸ)
        atmosphere_ocean = CrossRealmFluxes(momentum = momentum_flux_formulae)
    end

    return OceanSeaIceModelFluxes(surfaces,
                                  radiation,
                                  atmosphere_ocean,
                                  atmosphere_sea_ice,
                                  sea_ice_ocean)
end

Base.summary(crf::OceanSeaIceModelFluxes) = "OceanSeaIceModelFluxes"
Base.show(io::IO, crf::OceanSeaIceModelFluxes) = print(io, summary(crf))

#####
##### CrossRealmFlux
#####

struct RelativeAtmosphereOceanVelocity end
struct AtmosphereVelocity end

#####
##### Bulk formula
#####

struct BulkFormula{T, CD}
    transfer_velocity :: T
    transfer_coefficient :: CD
end

function BulkFormula(FT=Float64;
                     transfer_velocity = RelativeAtmosphereOceanVelocity(),
                     transfer_coefficient = 1e-3)

    return BulkFormula(transfer_velocity,
                       convert(FT, transfer_coefficient))
end

@inline Δϕt²(i, j, k, grid, ϕ1t, ϕ2, time) = @inbounds (ϕ1t[i, j, k, time] - ϕ2[i, j, k])^2

@inline function transfer_velocityᶠᶜᶜ(i, j, grid, time, ::RelativeAtmosphereOceanVelocity, Uₐ, Uₒ)
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    Δu = @inbounds uₐ[i, j, 1, time] - uₒ[i, j, 1]
    Δv² = ℑyᵃᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)
    return sqrt(Δu^2 + Δv²)
end

@inline function transfer_velocityᶜᶠᶜ(i, j, grid, time, ::RelativeAtmosphereOceanVelocity, Uₐ, Uₒ)
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v
    Δu² = ℑxᶜᵃᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv = @inbounds vₐ[i, j, 1, time] - vₒ[i, j, 1]
    return sqrt(Δu² + Δv^2)
end
