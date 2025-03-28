using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature
using ClimaSeaIce.SeaIceMomentumEquations: x_momentum_stress, y_momentum_stress

function compute_sea_ice_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice

    sea_ice_ocean_fluxes = coupled_model.interfaces.sea_ice_ocean_interface.fluxes
    interface_properties = coupled_model.interfaces.sea_ice_ocean_interface.properties 
   
    clock = sea_ice.model.clock

    Δt = ocean.Δt
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    Sᵢ = sea_ice.model.tracers.S
    ℵᵢ = sea_ice.model.ice_concentration
    hᵢ = sea_ice.model.ice_thickness
    h⁻ = coupled_model.interfaces.sea_ice_ocean_interface.previous_ice_thickness
    
    uᵢ, vᵢ = sea_ice.model.velocities

    dynamics = sea_ice.model.dynamics

    τs = if isnothing(dynamics)
        nothing
    else
        dynamics.external_momentum_stresses.bottom
    end

    ocean_properties = coupled_model.interfaces.ocean_properties
    liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)

    # What about the latent heat removed from the ocean when ice forms?
    # Is it immediately removed from the ocean? Or is it stored in the ice?
    launch!(arch, grid, :xy, _compute_sea_ice_ocean_fluxes!,
            sea_ice_ocean_fluxes, grid, clock, hᵢ, h⁻, ℵᵢ, Sᵢ, Tₒ, Sₒ, uᵢ, vᵢ, 
            τs, liquidus, ocean_properties, interface_properties, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_fluxes!(sea_ice_ocean_fluxes,
                                                grid,
                                                clock, 
                                                ice_thickness,
                                                previous_ice_thickness,
                                                ice_concentration,
                                                ice_salinity,
                                                ocean_temperature,
                                                ocean_salinity,
                                                sea_ice_u_velocity, 
                                                sea_ice_v_velocity, 
                                                sea_ice_ocean_stresses,
                                                liquidus,
                                                ocean_properties,
                                                interface_properties,
                                                Δt)

    i, j = @index(Global, NTuple)

    Nz  = size(grid, 3)
    Qᶠₒ = sea_ice_ocean_fluxes.frazil_heat
    Qᵢₒ = sea_ice_ocean_fluxes.interface_heat
    Jˢ  = sea_ice_ocean_fluxes.salt
    τx  = sea_ice_ocean_fluxes.x_momentum
    τy  = sea_ice_ocean_fluxes.y_momentum
    uᵢ  = sea_ice_u_velocity
    vᵢ  = sea_ice_v_velocity
    Tₒ  = ocean_temperature
    Sₒ  = ocean_salinity
    Sᵢ  = ice_salinity
    hᵢ  = ice_thickness
    ℵᵢ  = ice_concentration
    h⁻  = previous_ice_thickness
    ρₒ  = ocean_properties.reference_density
    cₒ  = ocean_properties.heat_capacity
    uₘ★ = interface_properties.characteristic_melting_speed

    δQ_frazil = zero(grid)

    for k = Nz:-1:1
        @inbounds begin
            # Various quantities
            Δz = Δzᶜᶜᶜ(i, j, k, grid)
            Tᵏ = Tₒ[i, j, k]
            Sᵏ = Sₒ[i, j, k]
        end

        # Melting / freezing temperature at the surface of the ocean
        Tₘ = melting_temperature(liquidus, Sᵏ)
        freezing = Tᵏ < Tₘ 

        # Compute change in ocean heat energy due to freezing.
        # When Tᵏ < Tₘ, we heat the ocean back to melting temperature by extracting heat from the ice.
        # The energy used to heat the ocean is transported instantneously to the surface and captured
        # by the sea ice componet (in reality, the formation of nascent ice crystals called frazil ice
        # transfers heat into the ocean).
        #
        δE_frazil = freezing * ρₒ * cₒ * (Tₘ - Tᵏ)

        # Perform temperature adjustment
        @inbounds Tₒ[i, j, k] = ifelse(freezing, Tₘ, Tᵏ)

        # Compute the heat flux from ocean into ice during frazil formation.
        #
        # A negative value δQ_frazil < 0 implies that heat is fluxed from the ice into
        # the ocean, cooling the ice and heating the ocean (δEₒ > 0). This occurs when
        # frazil ice is formed within the ocean.
        δQ_frazil -= δE_frazil * Δz / Δt
    end

    @inbounds begin
        Tᴺ = Tₒ[i, j, Nz]
        Sᴺ = Sₒ[i, j, Nz]
        ℵ  = ℵᵢ[i, j, 1]
    end

    # Compute total heat associated with temperature adjustment
    Tₘ = melting_temperature(liquidus, Sᴺ)
    δE_ice_bath = ρₒ * cₒ *  (Tₘ - Tᴺ)

    # Ice bath approximation
    adjust_temperature = (Tᴺ > Tₘ) & (ℵ > 0)
    @inbounds Tₒ[i, j, Nz] = ifelse(adjust_temperature, Tₘ, Tᴺ)

    # Compute the heat flux from ocean into ice due to sea ice melting.
    # A positive value δQ_melting > 0 corresponds to ocean cooling; ie
    # is fluxing upwards, into the ice. This occurs when applying the
    # ice bath equilibrium condition to cool down a warm ocean (δEₒ < 0).
    Δz = Δzᶜᶜᶜ(i, j, Nz, grid)
    δQ_melting = - δE_ice_bath * Δz / Δt #  uₘ★

    # Store column-integrated ice-ocean heat flux
    @inbounds Qᶠₒ[i, j, 1] = δQ_frazil
    @inbounds Qᵢₒ[i, j, 1] = δQ_melting * ℵ # Melting depends on concentration

    sea_ice_fields = (; u = uᵢ, v = vᵢ, h = hᵢ, ℵ = ℵᵢ)
    τₒᵢ = sea_ice_ocean_stresses

    @inbounds begin
        # Change in thickness
        Δh = hᵢ[i, j, 1] - h⁻[i, j, 1]

        # Update surface salinity flux.
        # Note: the Δt below is the ocean time-step, eg.
        # ΔS = ⋯ - ∮ Jˢ dt ≈ ⋯ - Δtₒ * Jˢ 
        Jˢ[i, j, 1] = Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz])

        # Update previous ice thickness
        h⁻[i, j, 1] = hᵢ[i, j, 1]
        τx[i, j, 1] = ℑyᵃᶜᵃ(i, j, 1, grid, x_momentum_stress, τₒᵢ, clock, sea_ice_fields)
        τy[i, j, 1] = ℑxᶜᵃᵃ(i, j, 1, grid, y_momentum_stress, τₒᵢ, clock, sea_ice_fields)
    end
end
