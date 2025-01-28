using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

function compute_sea_ice_ocean_fluxes!(coupled_model)
    #compute_sea_ice_ocean_salinity_flux!(coupled_model)
    compute_sea_ice_ocean_latent_heat_flux!(coupled_model)
    return nothing
end

function compute_sea_ice_ocean_latent_heat_flux!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    Qₒ = coupled_model.interfaces.sea_ice_ocean_interface.fluxes.heat
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
            Qₒ, grid, ℵᵢ, Tₒ, Sₒ, liquidus, ocean_properties, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_latent_heat_flux!(heat_flux,
                                                          grid,
                                                          ice_concentration,
                                                          ocean_temperature,
                                                          ocean_salinity,
                                                          liquidus,
                                                          ocean_properties,
                                                          Δt)

    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)
    Qᵢₒ = heat_flux
    Tₒ = ocean_temperature
    Sₒ = ocean_salinity
    ρₒ = ocean_properties.reference_density
    cₒ = ocean_properties.heat_capacity

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

    # Perform temperature adjustment at due to presence of sea ice
    kᴺ = size(grid, 3)
    Δz = Δzᶜᶜᶜ(i, j, kᴺ, grid)

    @inbounds begin
        Tᴺ = Tₒ[i, j, kᴺ]
        Sᴺ = Sₒ[i, j, kᴺ]
    end

    # Adjust temperature 
    Tₘ = melting_temperature(liquidus, Sᴺ)
    ΔT = ℵ * (Tₘ - Tᴺ)
    max_δQ = 1000
    max_δE = max_δQ * Δt / Δz
    max_ΔT = max_δE / (ρₒ * cₒ)
    ΔT = min(ΔT, + max_ΔT)
    ΔT = max(ΔT, - max_ΔT)
    @inbounds Tₒ[i, j, kᴺ] = Tᴺ + ΔT

    # Compute total heat associated with temperature adjustment
    δE_ice_bath = ρₒ * cₒ * ΔT

    # Compute the heat flux from ocean into ice due to sea ice melting.
    # A positive value δQ_melting > 0 corresponds to ocean cooling; ie
    # is fluxing upwards, into the ice. This occurs when applying the
    # ice bath equilibrium condition to cool down a warm ocean (δEₒ < 0).
    δQ_melting = - δE_ice_bath * Δz / Δt

    # @printf("Q_frazil: %.1f, Q_io: %.1f \n", δQ_frazil, δQ_melting)

    # Store column-integrated ice-ocean heat flux
    @inbounds Qᵢₒ[i, j, 1] = 0 #δQ_frazil + δQ_melting
end

function compute_sea_ice_ocean_salinity_flux!(coupled_model)
    # Compute salinity increment due to changes in ice thickness

    sea_ice = coupled_model.sea_ice
    ocean = coupled_model.ocean
    grid = ocean.model.grid
    arch = architecture(grid)
    Sₒ = ocean.model.tracers.S
    Sᵢ = sea_ice.model.ice_salinity
    Δt = ocean.Δt
    hⁿ = sea_ice.model.ice_thickness
    h⁻ = coupled_model.interfaces.previous_ice_thickness

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
    Jˢ = sea_ice_ocean_salinity_flux
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


