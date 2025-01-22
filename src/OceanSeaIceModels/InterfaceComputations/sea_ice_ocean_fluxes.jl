using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

function compute_sea_ice_ocean_fluxes!(coupled_model)
    #compute_sea_ice_ocean_salinity_flux!(coupled_model)
    compute_sea_ice_ocean_latent_heat_flux!(coupled_model)
    return nothing
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

function compute_sea_ice_ocean_latent_heat_flux!(coupled_model)
    ocean = coupled_model.ocean
    sea_ice = coupled_model.sea_ice
    ρₒ = coupled_model.interfaces.ocean_reference_density
    cₒ = coupled_model.interfaces.ocean_heat_capacity
    Qₒ = coupled_model.interfaces.sea_ice_ocean_interface.fluxes.heat
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    Δt = ocean.Δt
    hᵢ = sea_ice.model.ice_thickness

    liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)

    # What about the latent heat removed from the ocean when ice forms?
    # Is it immediately removed from the ocean? Or is it stored in the ice?
    launch!(arch, grid, :xy, _compute_sea_ice_ocean_latent_heat_flux!,
            Qₒ, grid, hᵢ, Tₒ, Sₒ, liquidus, ρₒ, cₒ, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_latent_heat_flux!(latent_heat_flux,
                                                          grid,
                                                          ice_thickness,
                                                          ocean_temperature,
                                                          ocean_salinity,
                                                          liquidus,
                                                          ρₒ, cₒ, Δt)

    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)
    Qₒ = latent_heat_flux
    hᵢ = ice_thickness
    Tₒ = ocean_temperature
    Sₒ = ocean_salinity

    δQ = zero(grid)
    icy_cell = @inbounds hᵢ[i, j, 1] > 0 # make ice bath approximation then

    for k = Nz:-1:1
        @inbounds begin
            # Various quantities
            Δz = Δzᶜᶜᶜ(i, j, k, grid)
            Tᴺ = Tₒ[i, j, k]
            Sᴺ = Sₒ[i, j, k]
        end

        # Melting / freezing temperature at the surface of the ocean
        Tₘ = melting_temperature(liquidus, Sᴺ)

        # Conditions for non-zero ice-ocean flux:
        #   - the ocean is below the freezing temperature, causing formation of ice.
        freezing = Tᴺ < Tₘ 

        #   - We are at the surface and the cell is covered by ice.
        icy_surface_cell = (k == Nz) & icy_cell

        # When there is a non-zero ice-ocean flux, we will instantaneously adjust the
        # temperature of the grid cells accordingly.
        adjust_temperature = freezing | icy_surface_cell

        # Compute change in ocean heat energy.
        #
        #   - When Tᴺ < Tₘ, we heat the ocean back to melting temperature by extracting heat from the ice,
        #     assuming that the heat flux (which is carried by nascent ice crystals called frazil ice) floats
        #     instantaneously to the surface.
        #
        #   - When Tᴺ > Tₘ and we are in a surface cell covered by ice, we assume equilibrium
        #     and cool the ocean by injecting excess heat into the ice.
        #
        δEₒ = adjust_temperature * ρₒ * cₒ * (Tₘ - Tᴺ)

        # Perform temperature adjustment
        @inline Tₒ[i, j, k] = ifelse(adjust_temperature, Tₘ, Tᴺ)

        # Compute the heat flux from ocean into ice.
        #
        # A positive value δQ > 0 implies that the ocean is cooled; ie heat
        # is fluxing upwards, into the ice. This occurs when applying the
        # ice bath equilibrium condition to cool down a warm ocean (δEₒ < 0).
        #
        # A negative value δQ < 0 implies that heat is fluxed from the ice into
        # the ocean, cooling the ice and heating the ocean (δEₒ > 0). This occurs when
        # frazil ice is formed within the ocean.

        δQ -= δEₒ * Δz / Δt
    end

    # Store ice-ocean flux
    @inbounds Qₒ[i, j, 1] = δQ
end
