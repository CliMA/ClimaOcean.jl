using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

function compute_sea_ice_ocean_fluxes!(coupled_model)
    compute_sea_ice_ocean_salinity_flux!(coupled_model)
    sea_ice_ocean_latent_heat_flux!(coupled_model)
    return nothing
end

function compute_sea_ice_ocean_salinity_flux!(coupled_model)
    # Compute salinity increment due to changes in ice thickness

    sea_ice = coupled_model.sea_ice
    ocean = coupled_model.ocean
    grid = ocean.model.grid
    arch = architecture(grid)
    Qˢ = ocean.model.tracers.S.boundary_conditions.top.condition
    Sₒ = ocean.model.tracers.S
    Sᵢ = sea_ice.model.tracers.S
    Δt = ocean.Δt
    hⁿ = sea_ice.model.ice_thickness
    h⁻ = coupled_model.fluxes.previous_ice_thickness
    ℵ  = sea_ice.model.ice_concentration

    liquidus = sea_ice.model.ice_thermodynamics.phase_transitions.liquidus
    ρₒ = coupled_model.fluxes.ocean_reference_density
    cₒ = coupled_model.fluxes.ocean_heat_capacity
    Qₒ = sea_ice.model.external_heat_fluxes.bottom
    Tₒ = ocean.model.tracers.T
    Δt = ocean.Δt

    # What about the latent heat removed from the ocean when ice forms?
    # Is it immediately removed from the ocean? Or is it stored in the ice?
    launch!(arch, grid, :xy, _compute_sea_ice_ocean_salinity_flux!,
            Qˢ, grid, ℵ, hⁿ, h⁻, Sᵢ, Sₒ, Δt)

    launch!(arch, grid, :xyz, _compute_sea_ice_ocean_latent_heat_flux!,
            Qₒ, grid, ℵ, hⁿ, Tₒ, Sₒ, ρₒ, cₒ, liquidus, Δt)

    return nothing
end

@kernel function _compute_sea_ice_ocean_salinity_flux!(sea_ice_ocean_salinity_flux,
                                                       grid,
                                                       ice_concentration,
                                                       ice_thickness,
                                                       previous_ice_thickness,
                                                       ice_salinity,
                                                       ocean_salinity,
                                                       Δt)
    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)

    hⁿ = ice_thickness
    h⁻ = previous_ice_thickness
    Qˢ = sea_ice_ocean_salinity_flux
    Sᵢ = ice_salinity
    Sₒ = ocean_salinity
    ℵ  = ice_concentration

    @inbounds begin
        # Change in thickness
        Δh = hⁿ[i, j, 1] - h⁻[i, j, 1]

        # Update surface salinity flux.
        # Note: the Δt below is the ocean time-step, eg.
        # ΔS = ⋯ - ∮ Qˢ dt ≈ ⋯ - Δtₒ * Qˢ 
        Qˢ[i, j, 1] = ifelse(ℵ > 0, Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz]), Qˢ[i, j, 1])

        # Update previous ice thickness
        h⁻[i, j, 1] = hⁿ[i, j, 1]
    end
end

@kernel function _compute_sea_ice_ocean_latent_heat_flux!(latent_heat_flux,
                                                          grid,
                                                          ice_concentration,
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
    ℵ  = ice_concentration

    δQ = zero(grid)
    icy_cell = @inbounds ℵ[i, j, 1] > 0 # make ice bath approximation then

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

function limit_fluxes_over_sea_ice!(grid, kernel_parameters, sea_ice::SeaIceModel,
                                    staggered_velocity_fluxes,
                                    net_tracer_fluxes,
                                    ocean_state)

    launch!(architecture(grid), grid, kernel_parameters, _massage_fluxes_on_sea_ice!,
            staggered_velocity_fluxes,
            net_tracer_fluxes,
            grid, 
            sea_ice.ice_concentration,
            sea_ice.ice_dynamics.ocean_ice_drag_coefficient,
            sea_ice.velocities,
            ocean_state.u,
            ocean_state.v)

    return nothing
end

@kernel function _massage_fluxes_on_sea_ice!(staggered_velocity_fluxes,
                                             net_tracer_fluxes,
                                             grid, ℵ, Cᴰ, 𝒰ᵢ, uₒ, vₒ)    

    i, j = @index(Global, NTuple)

    @inbounds begin
        Jᵘ = staggered_velocity_fluxes.u
        Jᵛ = staggered_velocity_fluxes.v
        Jᵀ = net_tracer_fluxes.T
        Jˢ = net_tracer_fluxes.S

        icy_cell = ℵ[i, j, 1] > 0
        cooling_sea_ice = icy_cell & (Jᵀ[i, j, 1] > 0)

        # If we are in a "sea ice" region we remove all fluxes
        Jᵘ[i, j, 1] = ifelse(icy_cell, Cᴰ * (uₒ[i, j, 1] - 𝒰ᵢ.u[i, j, 1]), Jᵘ[i, j, 1]) 
        Jᵛ[i, j, 1] = ifelse(icy_cell, Cᴰ * (vₒ[i, j, 1] - 𝒰ᵢ.v[i, j, 1]), Jᵛ[i, j, 1]) 
    end
end
