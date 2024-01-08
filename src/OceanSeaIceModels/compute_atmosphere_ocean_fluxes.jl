using Oceananigans.Grids: inactive_node
using Oceananigans.Fields: ConstantField

#####
##### Utilities
#####

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

#####
##### Computation
#####

const c = Center()
const f = Face()

function compute_atmosphere_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    atmosphere = coupled_model.atmosphere

    # Basic model properties
    grid = ocean.model.grid
    arch = architecture(grid)
    clock = ocean.model.clock

    # Ocean, atmosphere, and sea ice state
    ocean_velocities = surface_velocities(ocean)
    ocean_tracers    = surface_tracers(ocean)

    atmosphere_velocities            = surface_velocities(atmosphere)
    atmosphere_tracers               = surface_tracers(atmosphere)
    atmosphere_downwelling_radiation = downwelling_radiation(atmosphere)
    atmosphere_freshwater_flux       = freshwater_flux(atmosphere)

    ice_thickness = sea_ice_thickness(sea_ice)

    # Fluxes, and flux contributors
    net_momentum_fluxes         = coupled_model.surfaces.ocean.momentum
    net_tracer_fluxes           = coupled_model.surfaces.ocean.tracers
    bulk_momentum_flux_formulae = coupled_model.fluxes.atmosphere_ocean.momentum
    bulk_heat_flux_formulae     = coupled_model.fluxes.atmosphere_ocean.heat
    bulk_tracer_flux_formulae   = coupled_model.fluxes.atmosphere_ocean.tracers
    bulk_velocity_scale         = coupled_model.fluxes.bulk_velocity_scale
    surface_radiation           = coupled_model.fluxes.surface_radiation

    ocean_state = merge(ocean_velocities, ocean_tracers)

    atmosphere_state = (ρ = density(atmosphere),
                        cₚ = specific_heat(atmosphere))

    atmosphere_state = merge(atmosphere_velocities,
                             atmosphere_tracers,
                             atmosphere_state)

    launch!(arch, grid, :xy, _compute_atmosphere_ocean_fluxes!,
            grid, clock,
            net_momentum_fluxes,
            net_tracer_fluxes,
            bulk_velocity_scale,
            bulk_momentum_flux_formulae,
            bulk_heat_flux_formulae,
            bulk_tracer_flux_formulae,
            ocean_state,
            atmosphere_state,
            atmosphere_downwelling_radiation,
            surface_radiation,
            atmosphere_freshwater_flux,
            coupled_model.ocean_reference_density,
            coupled_model.ocean_heat_capacity,
            ice_thickness)

    return nothing
end

@kernel function _compute_atmosphere_ocean_fluxes!(grid,
                                                   clock,
                                                   net_momentum_fluxes,
                                                   net_tracer_fluxes,
                                                   bulk_velocity,
                                                   bulk_momentum_flux_formulae,
                                                   bulk_heat_flux_formulae,
                                                   bulk_tracer_flux_formulae,
                                                   ocean_state,
                                                   atmos_state,
                                                   downwelling_radiation,
                                                   surface_radiation,
                                                   prescribed_freshwater_flux,
                                                   ocean_reference_density,
                                                   ocean_heat_capacity,
                                                   ice_thickness)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) # index of the top ocean cell

    time = Time(clock.time)

    Jᵘ = net_momentum_fluxes.u
    Jᵛ = net_momentum_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    # Note: there could one or more formula(e)
    τˣ_formula = bulk_momentum_flux_formulae.u
    τʸ_formula = bulk_momentum_flux_formulae.v
    Q_formula = bulk_heat_flux_formulae
    F_formula = bulk_tracer_flux_formulae.S

    atmos_state_names = keys(atmos_state)
    ocean_state_names = keys(atmos_state)

    atmos_state_ij = stateindex(atmos_state, i, j, 1, time)
    ocean_state_ij = stateindex(ocean_state, i, j, 1, time)

    # Compute transfer velocity scale
    ΔUᶠᶜᶜ = bulk_velocity_scaleᶠᶜᶜ(i, j, grid, time, bulk_velocity, atmos_state, ocean_state)
    ΔUᶜᶠᶜ = bulk_velocity_scaleᶜᶠᶜ(i, j, grid, time, bulk_velocity, atmos_state, ocean_state)
    ΔUᶜᶜᶜ = bulk_velocity_scaleᶜᶜᶜ(i, j, grid, time, bulk_velocity, atmos_state, ocean_state)

    # Compute momentum fluxes
    τˣ = cross_realm_flux(i, j, grid, time, τˣ_formula, ΔUᶠᶜᶜ, atmos_state, ocean_state)
    τʸ = cross_realm_flux(i, j, grid, time, τʸ_formula, ΔUᶜᶠᶜ, atmos_state, ocean_state)

    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, downwelling_radiation, surface_radiation)
     
    Qu = net_upwelling_radiation(i, j, grid, time, surface_radiation, ocean_state)
    Q★ = cross_realm_flux(i, j, grid, time, Q_formula, ΔUᶜᶜᶜ, atmos_state_ij, ocean_state_ij)
    Q = Q★ + Qd + Qu

    # Compute salinity fluxes, bulk flux first
    Fp = cross_realm_flux(i, j, grid, time, prescribed_freshwater_flux)
    F★ = cross_realm_flux(i, j, grid, time, F_formula, ΔUᶜᶜᶜ, atmos_state_ij, ocean_state_ij)
    F = F★ + Fp

    # Then the rest of the heat fluxes
    ρₒ = ocean_reference_density
    cₚ = ocean_heat_capacity

    atmos_ocean_Jᵘ = τˣ / ρₒ
    atmos_ocean_Jᵛ = τʸ / ρₒ
    atmos_ocean_Jᵀ = Q / (ρₒ * cₚ)

    S = ocean_state_ij.S
    atmos_ocean_Jˢ = S * F

    @inbounds begin
        # Set fluxes
        # TODO: should this be peripheral_node?
        Jᵘ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, f, c, c), zero(grid), atmos_ocean_Jᵘ)
        Jᵛ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, f, c), zero(grid), atmos_ocean_Jᵛ)
        Jᵀ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, c, c), zero(grid), atmos_ocean_Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive_node(i, j, kᴺ, grid, c, c, c), zero(grid), atmos_ocean_Jˢ)
    end
end

@inline function net_downwelling_radiation(i, j, grid, time, downwelling_radiation, surface_radiation)
    Qˢʷ = downwelling_radiation.shortwave
    Qˡʷ = downwelling_radiation.longwave
    α = stateindex(surface_radiation.reflection.ocean, i, j, 1, time)

    return @inbounds - (1 - α) * Qˢʷ[i, j, 1, time] - Qˡʷ[i, j, 1, time]
end

@inline function net_upwelling_radiation(i, j, grid, time, surface_radiation, ocean_state)
    σ = surface_radiation.stefan_boltzmann_constant
    Tᵣ = surface_radiation.reference_temperature
    ϵ = stateindex(surface_radiation.emission.ocean, i, j, 1, time)

    # Ocean surface temperature (departure from reference, typically in ᵒC)
    Tₒ = @inbounds ocean_state.T[i, j, 1]

    # Note: positive implies _upward_ heat flux, and therefore cooling.
    return σ * ϵ * (Tₒ + Tᵣ)^4
end

