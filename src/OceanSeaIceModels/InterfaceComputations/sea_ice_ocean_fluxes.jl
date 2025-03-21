using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

function compute_sea_ice_ocean_fluxes!(coupled_model)
    compute_sea_ice_ocean_salinity_flux!(coupled_model)
    compute_sea_ice_ocean_latent_heat_flux!(coupled_model)
    return nothing
end

function compute_sea_ice_ocean_latent_heat_flux!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    Qᶠₒ = coupled_model.interfaces.sea_ice_ocean_interface.fluxes.frazil_heat
    Qᵢₒ = coupled_model.interfaces.sea_ice_ocean_interface.fluxes.interface_heat
    
    interface_properties = coupled_model.interfaces.sea_ice_ocean_interface.properties 
   
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    Δt = ocean.Δt
    ℵᵢ = sea_ice.model.ice_concentration
    
    ocean_properties = coupled_model.interfaces.ocean_properties
    liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)

    # What about the latent heat removed from the ocean when ice forms?
    # Is it immediately removed from the ocean? Or is it stored in the ice?
    launch!(arch, grid, :xy, _compute_sea_ice_ocean_latent_heat_flux!,
            Qᶠₒ, Qᵢₒ, grid, ℵᵢ, Tₒ, Sₒ, liquidus, ocean_properties, interface_properties, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_latent_heat_flux!(frazil_heat_flux,
                                                          interface_heat_flux,
                                                          grid,
                                                          ice_concentration,
                                                          ocean_temperature,
                                                          ocean_salinity,
                                                          liquidus,
                                                          ocean_properties,
                                                          interface_properties,
                                                          Δt)

    i, j = @index(Global, NTuple)

    Nz  = size(grid, 3)
    Qᶠₒ = frazil_heat_flux
    Qᵢₒ = interface_heat_flux
    Tₒ  = ocean_temperature
    Sₒ  = ocean_salinity
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

    # Ice bath approximation
    adjust_temperature = (Tᴺ > Tₘ) & (ℵ > 0)
    @inbounds Tₒ[i, j, Nz] = ifelse(adjust_temperature, Tₘ, Tᴺ)

    # Compute the heat flux from ocean into ice due to sea ice melting.
    # A positive value δQ_melting > 0 corresponds to ocean cooling; ie
    # is fluxing upwards, into the ice. This occurs when applying the
    # ice bath equilibrium condition to cool down a warm ocean (δEₒ < 0).
    Δz = Δzᶜᶜᶜ(i, j, Nz, grid)
    δQ_melting = - δE_ice_bath * Δz / Δt  

    # Store column-integrated ice-ocean heat flux
    @inbounds Qᶠₒ[i, j, 1] = δQ_frazil
    @inbounds Qᵢₒ[i, j, 1] = δQ_melting * ℵ # Melting depends on concentration
end

function compute_sea_ice_ocean_salinity_flux!(coupled_model)
    # Compute salinity increment due to changes in ice thickness

    sea_ice = coupled_model.sea_ice
    ocean = coupled_model.ocean
    grid = sea_ice.model.grid
    arch = architecture(grid)
    Sₒ = ocean.model.tracers.S
    Sᵢ = sea_ice.model.tracers.S
    Δt = ocean.Δt
    hⁿ = sea_ice.model.ice_thickness
    h⁻ = coupled_model.interfaces.sea_ice_ocean_interface.previous_ice_thickness

    interface_fluxes = coupled_model.interfaces.sea_ice_ocean_interface.fluxes

    launch!(arch, grid, :xy, _compute_sea_ice_ocean_salinity_flux!,
            interface_fluxes.salt, grid, hⁿ, h⁻, Sᵢ, Sₒ, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_salinity_flux!(salt_flux,
                                                       grid,
                                                       ice_thickness,
                                                       previous_ice_thickness,
                                                       ice_salinity,
                                                       ocean_salinity,
                                                       Δt)
    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)

    hⁿ = ice_thickness
    h⁻ = previous_ice_thickness
    Sᵢ = ice_salinity
    Sₒ = ocean_salinity
    Jˢ = salt_flux

    @inbounds begin
        # Change in thickness
        Δh = hⁿ[i, j, 1] - h⁻[i, j, 1]

        # Update surface salinity flux.
        # Note: the Δt below is the ocean time-step, eg.
        # ΔS = ⋯ - ∮ Jˢ dt ≈ ⋯ - Δtₒ * Jˢ 
        Jˢ[i, j, 1] = Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz])

        # Update previous ice thickness
        h⁻[i, j, 1] = hⁿ[i, j, 1]
    end
end
