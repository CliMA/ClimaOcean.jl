using Oceananigans.Grids: inactive_node

#####
##### Utilities
#####

# surface_velocities(ocean::Simulation{<:HydrostaticFreeSurfaceModel}) = _surface_velocities(ocean)

function surface_velocities(ocean::Simulation{<:HydrostaticFreeSurfaceModel})
    grid = ocean.model.grid
    Nz = size(grid, 3)
    u = interior(ocean.model.velocities.u, :, :, Nz)
    v = interior(ocean.model.velocities.v, :, :, Nz)
    w = interior(ocean.model.velocities.w, :, :, Nz+1)
    return (; u, v, w)
end

function surface_tracers(ocean::Simulation{<:HydrostaticFreeSurfaceModel})
    grid = ocean.model.grid
    Nz = size(grid, 3)
    tracers = ocean.model.tracers
    tracer_names = keys(tracers)
    sfc_tracers = NamedTuple(name => interior(tracers[name], :, :, Nz)
                             for name in tracer_names)
    return sfc_tracers
end

#####
##### Computation
#####

const c = Center()
const f = Face()

function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    atmosphere = coupled_model.atmosphere

    # Basic model properties
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = ocean.model.clock

    # Ocean, atmosphere, and sea ice state
    ocean_velocities = surface_velocities(ocean)
    ocean_tracers    = surface_tracers(ocean)

    atmosphere_velocities = surface_velocities(atmosphere)
    atmosphere_tracers    = nothing

    ice_thickness = coupled_model.sea_ice.model.ice_thickness

    # Fluxes, and flux contributors
    net_momentum_fluxes         = coupled_model.fluxes.surfaces.ocean.momentum
    net_tracer_fluxes           = coupled_model.fluxes.surfaces.ocean.tracers
    radiation                   = coupled_model.fluxes.radiation
    momentum_flux_contributions = coupled_model.fluxes.atmosphere_ocean.momentum
    heat_flux_contributions     = coupled_model.fluxes.atmosphere_ocean.heat
    tracer_flux_contributions   = coupled_model.fluxes.atmosphere_ocean.heat

    # Parameters?
    atmosphere_ocean_parameters = (
        ρₐ = 1.2,
        ρₒ = coupled_model.ocean_reference_density,
        cₚ = 3991.0,
    )

    launch!(arch, grid, :xy, _compute_atmosphere_ocean_fluxes!,
            grid, clock,
            net_momentum_fluxes,
            net_tracer_fluxes,
            momentum_flux_contributions,
            heat_flux_contributions,
            tracer_flux_contributions,
            ocean_velocities,
            atmosphere_velocities,
            ocean_tracers,
            atmosphere_tracers,
            radiation,
            atmosphere_ocean_parameters,
            ice_thickness)

    return nothing
end

@inline function air_sea_difference(i, j, grid, time, ::RelativeVelocityScale,
                                    air::SomeKindOfFieldTimeSeries, sea::AbstractArray)

    δ = @inbounds air[i, j, 1, time] - sea[i, j, 1]
    return δ
end

@kernel function _compute_atmosphere_ocean_fluxes!(grid,
                                                   clock,
                                                   net_momentum_fluxes,
                                                   net_tracer_fluxes,
                                                   momentum_flux_contributions,
                                                   heat_flux_contributions,
                                                   tracer_flux_contributions,
                                                   ocean_velocities,
                                                   atmosphere_velocities,
                                                   ocean_tracers,
                                                   atmosphere_tracers,
                                                   radiation,
                                                   atmosphere_ocean_parameters,
                                                   ice_thickness)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) # index of the top ocean cell

    time = Time(clock.time)

    τˣ = net_momentum_fluxes.u
    τʸ = net_momentum_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    Uₒ = ocean_velocities
    Uₐ = atmosphere_velocities
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    Tₒ = ocean_tracers.T

    ρₐ = atmosphere_ocean_parameters.ρₐ
    ρₒ = atmosphere_ocean_parameters.ρₒ
    cₚ = atmosphere_ocean_parameters.cₚ

    u_formula = momentum_flux_contributions.u.velocity_scale
    v_formula = momentum_flux_contributions.v.velocity_scale
    cᴰ = momentum_flux_contributions.u.transfer_coefficient

    # Compute transfer velocity scale
    Vᶠᶜᶜ = bulk_velocity_scaleᶠᶜᶜ(i, j, grid, time, u_formula, Uₐ, Uₒ)
    Vᶜᶠᶜ = bulk_velocity_scaleᶜᶠᶜ(i, j, grid, time, v_formula, Uₐ, Uₒ)

    Δu = air_sea_difference(i, j, grid, time, u_formula, uₐ, uₒ)
    Δv = air_sea_difference(i, j, grid, time, v_formula, vₐ, vₒ)

    @inbounds begin
        atmos_ocean_τˣ = - ρₐ / ρₒ * cᴰ * Δu * Vᶠᶜᶜ
        atmos_ocean_τʸ = - ρₐ / ρₒ * cᴰ * Δv * Vᶜᶠᶜ

        # TODO: should this be peripheral_node?
        τˣ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, f, c, c), zero(grid), atmos_ocean_τˣ)
        τʸ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, f, c), zero(grid), atmos_ocean_τʸ)

        # Radiation first
        Q  = downwelling_radiation(i, j, grid, time, radiation)
        Q += upwelling_radiation(i, j, grid, time, radiation, Tₒ)

        # Then the rest of the heat fluxes
        atmos_ocean_Jᵀ = Q / (ρₒ * cₚ)
        Jᵀ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, c, c), zero(grid), atmos_ocean_Jᵀ)
    end
end

@inline Δϕt²(i, j, k, grid, ϕ1t, ϕ2, time) = @inbounds (ϕ1t[i, j, k, time] - ϕ2[i, j, k])^2

@inline function downwelling_radiation(i, j, grid, time, radiation)
    Qˢʷ = radiation.downwelling_shortwave_radiation
    Qˡʷ = radiation.downwelling_longwave_radiation
    α = radiation.ocean_albedo
    return @inbounds - (1 - α) * Qˢʷ[i, j, 1, time] - Qˡʷ[i, j, 1, time]
end

@inline function upwelling_radiation(i, j, grid, time, radiation, Tₒ)
    σ = radiation.stefan_boltzmann_constant
    ϵ = radiation.ocean_emissivity
    Tᵣ = radiation.reference_temperature

    # Note: positive implies _upward_ heat flux, and therefore cooling.
    return @inbounds σ * ϵ * (Tₒ[i, j, 1] + Tᵣ)^4
end

@inline function bulk_velocity_scaleᶠᶜᶜ(i, j, grid, time, ::RelativeVelocityScale, Uₐ, Uₒ)
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    Δu = @inbounds uₐ[i, j, 1, time] - uₒ[i, j, 1]
    Δv² = ℑyᵃᶜᵃ(i, j, 1, grid, Δϕt², vₐ, vₒ, time)

    return sqrt(Δu^2 + Δv²)
end

@inline function bulk_velocity_scaleᶜᶠᶜ(i, j, grid, time, ::RelativeVelocityScale, Uₐ, Uₒ)
    uₐ = Uₐ.u
    vₐ = Uₐ.v
    uₒ = Uₒ.u
    vₒ = Uₒ.v

    Δu² = ℑxᶜᵃᵃ(i, j, 1, grid, Δϕt², uₐ, uₒ, time)
    Δv = @inbounds vₐ[i, j, 1, time] - vₒ[i, j, 1]

    return sqrt(Δu² + Δv^2)
end

