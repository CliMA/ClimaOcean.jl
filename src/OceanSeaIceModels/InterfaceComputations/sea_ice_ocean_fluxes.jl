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
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    ℵᵢ = sea_ice.model.ice_concentration
    hᵢ = sea_ice.model.ice_thickness

    Q = coupled_model.interfaces.sea_ice_ocean_interface.fluxes.heat
    ocean_properties = coupled_model.interfaces.ocean_properties
    sea_ice_properties = coupled_model.interfaces.sea_ice_properties
    phase_transitions = sea_ice.model.ice_thermodynamics.phase_transitions
    grid = ocean.model.grid
    arch = architecture(grid)

    Ch = coupled_model.interfaces.flux_formulation.heat_transfer_coefficient

    # What about the latent heat removed from the ocean when ice forms?
    # Is it immediately removed from the ocean? Or is it stored in the ice?
    launch!(arch, grid, :xy, _compute_sea_ice_ocean_latent_heat_flux!,
            Q, grid, hᵢ, ℵᵢ, Tₒ, Sₒ, Ch, phase_transitions, ocean_properties, sea_ice_properties)

    return nothing
end

@kernel function _compute_sea_ice_ocean_latent_heat_flux!(heat_flux,
                                                          grid,
                                                          ice_thickness,
                                                          ice_concentration,
                                                          ocean_temperature,
                                                          ocean_salinity,
                                                          heat_transfer_coefficient,
                                                          phase_transitions,
                                                          ocean_properties,
                                                          sea_ice_properties)

    i, j = @index(Global, NTuple)

    Nz  = size(grid, 3)
    Qᵢₒ = heat_flux
    Tₒ  = ocean_temperature
    Sₒ  = ocean_salinity
    ρₒ  = ocean_properties.reference_density
    cₒ  = ocean_properties.heat_capacity
    ρᵢ  = sea_ice_properties.reference_density
    cᵢ  = sea_ice_properties.heat_capacity
    Ch  = heat_transfer_coefficients
    liquidus = phase_transitions.liquidus

    ℵ = @inbounds ice_concentration[i, j, 1]
    h = @inbounds ice_thickness[i, j, 1]

    ΔVᵢ = zero(grid)

    for k = Nz:-1:1
        @inbounds begin
            # Various quantities
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
        # Perform temperature adjustment
        @inbounds Tₒ[i, j, k] = ifelse(freezing, Tₘ, Tᵏ)

        # Compute the change in ice mass due to frazil ice formation
        # following conservation of energy:
        # mₒ ΔE = mᵢ (Eᵢ - mᵢ Ef)
        ℰₘ = SeaIceThermodynamics.latent_heat(phase_transitions, Tₘ)
        
        ΔT = (Tᵏ - Tₘ) * freezing
        Vₒ = Vᶜᶜᶜ(i, j, k, grid) * ρₒ
        Tₘ = convert_to_kelvin(ocean_properties.temperature_units, Tₘ)

        # Change in (per volume) energy of sea ice
        ΔEᵢ = ρᵢ * Tᵐ * (cᵢ - cₒ) - ℰₘ
        
        # Change in volume of sea ice due to frazil ice formation
        ΔVᵢ += Vₒ * ρₒ * cₒ * ΔT / ΔEᵢ
    end

    # Compute the thickness change in the sea ice due to the frazil ice formation 
    @inbounds h[i, j, 1] += ΔVᵢ / Azᶜᶜᶜ(i, j, Nz, grid)

    @inbounds begin
        Tᴺ = Tₒ[i, j, Nz]
        Sᴺ = Sₒ[i, j, Nz]
    end

    Tₘ = melting_temperature(liquidus, Sᴺ)
    ΔT = ℵ * (Tₘ - Tᴺ)

    # Compute a heat flux, based on a simple turbulent heat transfer 
    @inbounds Qᵢₒ[i, j, 1] = ρₒ * cₒ * Ch * ΔT
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


