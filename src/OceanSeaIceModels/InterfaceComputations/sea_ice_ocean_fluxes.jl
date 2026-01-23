using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaOcean.OceanSeaIceModels: ocean_temperature, ocean_salinity
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using ClimaSeaIce.SeaIceDynamics: x_momentum_stress, y_momentum_stress

"""
    compute_sea_ice_ocean_fluxes!(coupled_model)

Compute heat, salt, and momentum fluxes at the sea ice-ocean interface.

This function computes:
- Frazil heat flux: heat released when ocean temperature drops below freezing (all formulations)
- Interface heat flux: heat flux from ocean to ice, computed using the specified formulation
- Salt flux: salt exchange due to ice growth/melt
- Momentum stresses: ice-ocean momentum transfer

The interface heat flux formulation is determined by `coupled_model.interfaces.sea_ice_ocean_interface.flux_formulation`.
"""
function compute_sea_ice_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    interface = coupled_model.interfaces.sea_ice_ocean_interface
    ocean_properties = coupled_model.interfaces.ocean_properties

    compute_sea_ice_ocean_fluxes!(interface, ocean, sea_ice, ocean_properties)

    return nothing
end

function compute_sea_ice_ocean_fluxes!(interface, ocean, sea_ice, ocean_properties)
    Δt = sea_ice.Δt
    Tₒ = ocean_temperature(ocean)
    Sₒ = ocean_salinity(ocean)
    Sⁱ = sea_ice.model.tracers.S
    ℵ = sea_ice.model.ice_concentration
    hᵢ = sea_ice.model.ice_thickness
    hc = sea_ice.model.ice_consolidation_thickness

    phase_transitions = sea_ice.model.ice_thermodynamics.phase_transitions
    liquidus = phase_transitions.liquidus
    L = phase_transitions.reference_latent_heat

    grid = sea_ice.model.grid
    clock = sea_ice.model.clock
    arch = architecture(grid)

    uᵢ, vᵢ = sea_ice.model.velocities
    dynamics = sea_ice.model.dynamics

    # Get interface data
    fluxes = interface.fluxes
    flux_formulation = interface.flux_formulation
    Tᵢ = interface.temperature
    Sᵢ = interface.salinity

    if !isnothing(dynamics)
        kernel_parameters = interface_kernel_parameters(grid)
        τₛ = dynamics.external_momentum_stresses.bottom
        launch!(arch, grid, kernel_parameters, _compute_sea_ice_ocean_stress!,
                fluxes, grid, clock, hᵢ, ℵ, uᵢ, vᵢ, τₛ)
    else
        τₛ = nothing
    end

    launch!(arch, grid, :xy, _compute_sea_ice_ocean_fluxes!,
            flux_formulation, fluxes, Tᵢ, Sᵢ, grid, clock,
            hᵢ, hc, ℵ, Sⁱ, Tₒ, Sₒ, uᵢ, vᵢ, τₛ,
            liquidus, ocean_properties, L, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_stress!(fluxes, 
                                                grid, 
                                                clock, 
                                                ice_thickness,
                                                ice_concentration,
                                                sea_ice_u_velocity,
                                                sea_ice_v_velocity,
                                                sea_ice_ocean_stress)
    i, j = @index(Global, NTuple)

    τˣ = fluxes.x_momentum
    τʸ = fluxes.y_momentum
    Nz = size(grid, 3)
    
    uᵢ = sea_ice_u_velocity
    vᵢ = sea_ice_v_velocity
    hᵢ = ice_thickness
    ℵ = ice_concentration
    sea_ice_fields = (; u = uᵢ, v = vᵢ, h = hᵢ, ℵ = ℵ)

    # Momentum stresses
    @inbounds begin
        τˣ[i, j, 1] = x_momentum_stress(i, j, Nz, grid, sea_ice_ocean_stress, clock, sea_ice_fields)
        τʸ[i, j, 1] = y_momentum_stress(i, j, Nz, grid, sea_ice_ocean_stress, clock, sea_ice_fields)
    end
end

@kernel function _compute_sea_ice_ocean_fluxes!(flux_formulation,
                                                fluxes,
                                                interface_temperature,
                                                interface_salinity,
                                                grid,
                                                clock,
                                                ice_thickness,
                                                ice_consolidation_thickness,
                                                ice_concentration,
                                                ice_salinity,
                                                ocean_temperature,
                                                ocean_salinity,
                                                sea_ice_u_velocity,
                                                sea_ice_v_velocity,
                                                sea_ice_ocean_stresses,
                                                liquidus,
                                                ocean_properties,
                                                latent_heat,
                                                Δt)

    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)
    Qᶠ = fluxes.frazil_heat
    Qᵢ = fluxes.interface_heat
    Jˢ = fluxes.salt
    τˣ = fluxes.x_momentum
    τʸ = fluxes.y_momentum
    T★ = interface_temperature
    S★ = interface_salinity
    Tₒ = ocean_temperature
    Sₒ = ocean_salinity
    hc = ice_consolidation_thickness
    ℰ  = latent_heat

    ρₒ = ocean_properties.reference_density
    cₒ = ocean_properties.heat_capacity

    # =============================================
    # Part 1: Frazil ice formation (all formulations)
    # =============================================
    # When ocean temperature drops below freezing, frazil ice forms
    # and heat is released to the ice component.

    δQᶠ = zero(grid)

    for k = Nz:-1:1
        @inbounds begin
            Δz = Δzᶜᶜᶜ(i, j, k, grid)
            Tᵏ = Tₒ[i, j, k]
            Sᵏ = Sₒ[i, j, k]
        end

        # Melting/freezing temperature at this depth
        Tₘ = melting_temperature(liquidus, Sᵏ)
        freezing = Tᵏ < Tₘ

        # Compute change in ocean heat energy due to freezing.
        # When Tᵏ < Tₘ, we heat the ocean back to melting temperature
        # by extracting heat from the ice.
        δE = freezing * ρₒ * cₒ * (Tₘ - Tᵏ)

        # Perform temperature adjustment
        @inbounds Tₒ[i, j, k] = ifelse(freezing, Tₘ, Tᵏ)

        # Compute the heat flux from ocean into ice during frazil formation.
        # A negative value δQᶠ < 0 implies heat is fluxed from the ice into
        # the ocean (frazil ice formation).
        δQᶠ -= δE * Δz / Δt
    end

    # Store frazil heat flux
    @inbounds Qᶠ[i, j, 1] = δQᶠ

    # Freezing rate
    qᶠ = δQᶠ / ℰ

    @inbounds begin
        Tᴺ = Tₒ[i, j, Nz]               
        Sᴺ = Sₒ[i, j, Nz]               
        Sᵢ = ice_salinity[i, j, 1]      
        hᵢ = ice_thickness[i, j, 1]     
        ℵᵢ = ice_concentration[i, j, 1] 
        hc = ice_consolidation_thickness[i, j, 1] 
    end

    # Extract internal temperature (for ConductiveFluxTEF, zero otherwise)
    Tᵢ = extract_internal_temperature(flux_formulation, i, j)

    # Package states
    ocean_surface_state = (; T = Tᴺ, S = Sᴺ)
    ice_state = (; S = Sᵢ, h = hᵢ, hc = hc, ℵ = ℵᵢ, T = Tᵢ)

    # Compute friction velocity
    u★ = get_friction_velocity(flux_formulation.friction_velocity, i, j, grid, τˣ, τʸ, ρₒ)

    # =============================================
    # Part 3: Interface heat flux (formulation-specific)
    # =============================================
    # Returns interfacial heat flux Q, melt rate qᵐ, and interface T, S
    Qᵢₒ, qᵐ, Tᵦ, Sᵦ = compute_interface_heat_flux(flux_formulation,
                                                   ocean_surface_state, ice_state,
                                                   liquidus, ocean_properties, ℰ, u★)

    # Store interface values and heat flux
    @inbounds T★[i, j, 1] = Tᵦ
    @inbounds S★[i, j, 1] = Sᵦ
    @inbounds Qᵢ[i, j, 1] = Qᵢₒ

    # =============================================
    # Part 4: Salt flux
    # =============================================
    # Salt flux from melting/freezing:
    # - during ice melt   (qᵐ > 0), fresh meltwater dilutes the ocean
    # - during ice growth (qᶠ < 0), brine rejection adds salt to ocean
    @inbounds Jˢ[i, j, 1] = (qᵐ + qᶠ) / ρₒ * (Sᴺ - Sᵢ)
end
