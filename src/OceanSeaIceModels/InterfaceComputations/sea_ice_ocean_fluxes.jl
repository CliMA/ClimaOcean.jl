using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

function compute_sea_ice_ocean_fluxes!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    fluxes = coupled_model.interfaces.sea_ice_ocean_interface.fluxes
    
    interface_properties = coupled_model.interfaces.sea_ice_ocean_interface.properties 
   
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    Sᵢ = sea_ice.model.tracers.S
    hᵢ = sea_ice.model.ice_thickness
    h⁻ = coupled_model.interfaces.sea_ice_ocean_interface.previous_ice_thickness
    Δt = ocean.Δt
    ℵᵢ = sea_ice.model.ice_concentration
    
    ocean_properties = coupled_model.interfaces.ocean_properties
    liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)

    # What about the latent heat removed from the ocean when ice forms?
    # Is it immediately removed from the ocean? Or is it stored in the ice?
    launch!(arch, grid, :xy, _compute_sea_ice_ocean_fluxes!,
    fluxes, grid, hᵢ, ℵᵢ, Sᵢ, h⁻, Tₒ, Sₒ, liquidus, ocean_properties, interface_properties, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_fluxes!(sea_ice_ocean_fluxes,
                                                grid,
                                                ice_thickness,
                                                ice_concentration,
                                                ice_salinity,
                                                previous_ice_thickness,
                                                ocean_temperature,
                                                ocean_salinity,
                                                liquidus,
                                                ocean_properties,
                                                interface_properties,
                                                Δt)

    i, j = @index(Global, NTuple)

    Nz  = size(grid, 3)
    Qᶠₒ = sea_ice_ocean_fluxes.frazil_heat
    Qᵢₒ = sea_ice_ocean_fluxes.interface_heat
    Jˢ  = sea_ice_ocean_fluxes.salt
    Tₒ  = ocean_temperature
    Sₒ  = ocean_salinity
    Sᵢ  = ice_salinity
    hⁿ  = ice_thickness
    h⁻  = previous_ice_thickness
    ρₒ  = ocean_properties.reference_density
    cₒ  = ocean_properties.heat_capacity
    uₘ★ = interface_properties.characteristic_melting_speed

    ℵ = @inbounds ice_concentration[i, j, 1]
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
    end

    # Compute total heat associated with temperature adjustment
    Tₘ = melting_temperature(liquidus, Sᴺ)
    δE_ice_bath = ρₒ * cₒ *  (Tₘ - Tᴺ)

    # Compute the heat flux from ocean into ice due to sea ice melting.
    # A positive value δQ_melting > 0 corresponds to ocean cooling; ie
    # is fluxing upwards, into the ice. This occurs when applying the
    # ice bath equilibrium condition to cool down a warm ocean (δEₒ < 0).
    δQ_melting = - δE_ice_bath * uₘ★  

    @inbounds begin
        # Change in thickness
        Δh = hⁿ[i, j, 1] - h⁻[i, j, 1]

        # Update previous ice thickness
        h⁻[i, j, 1] = hⁿ[i, j, 1]

        # Update surface salinity flux.
        # Note: the Δt below is the ocean time-step, eg.
        # ΔS = ⋯ - ∮ Jˢ dt ≈ ⋯ - Δtₒ * Jˢ 
        Jˢ[i, j, 1] = Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz])

        # Update heat fluxes.
        Qᶠₒ[i, j, 1] = δQ_frazil
        Qᵢₒ[i, j, 1] = δQ_melting * ℵ # Melting depends on concentration
    end
end