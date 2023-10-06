
# If there is no atmosphere, do not compute fluxes! (this model has the ocean component which 
# will handle the top boundary_conditions, for example if we want to impose a value BC)
compute_air_sea_flux!(coupled_model::NoAtmosphereModel) = nothing

function compute_air_sea_flux!(coupled_model)
    ocean   = coupled_model.ocean
    forcing = coupled_model.atmospheric_forcing

    (; T, S) = ocean.model.tracers
    (; u, v) = ocean.model.velocities

    grid   = ocean.model.grid
    clock  = ocean.model.clock
    fields = prognostic_fields(ocean.model)

    Qˢ = T.boundary_conditions.top.condition
    Fˢ = S.boundary_conditions.top.condition
    τˣ = u.boundary_conditions.top.condition
    τʸ = v.boundary_conditions.top.condition

    ε  = coupled_model.ocean_emissivity
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    I₀ = coupled_model.solar_insolation

    ice_thickness = coupled_model.ice.model.ice_thickness

    launch!(ocean, :xy, _calculate_air_sea_fluxes!, Qˢ, Fˢ, τˣ, τʸ, ρₒ, cₒ, ε, I₀, 
            grid, clock, fields, forcing, ice_thickness)

    return nothing
end

function compute_ice_ocean_flux!(coupled_model)

    # probably need to expand this
    compute_ice_ocean_salinity_flux!(coupled_model)
    ice_ocean_latent_heat!(coupled_model)

    return nothing
end

function compute_ice_ocean_salinity_flux!(coupled_model)
    # Compute salinity increment due to changes in ice thickness

    ice = coupled_model.ice
    ocean = coupled_model.ocean
    grid = ocean.model.grid
    arch = architecture(grid)
    Qˢ = ocean.model.tracers.S.boundary_conditions.top.condition
    Sₒ = ocean.model.tracers.S
    Sᵢ = ice.model.ice_salinity
    Δt = ocean.Δt
    hⁿ = ice.model.ice_thickness
    h⁻ = coupled_model.previous_ice_thickness

    launch!(arch, grid, :xy, _compute_ice_ocean_salinity_flux!,
            Qˢ, grid, hⁿ, h⁻, Sᵢ, Sₒ, Δt)

    return nothing
end

@kernel function _compute_ice_ocean_salinity_flux!(ice_ocean_salinity_flux,
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
    Qˢ = ice_ocean_salinity_flux
    Sᵢ = ice_salinity
    Sₒ = ocean_salinity

    @inbounds begin
        # Thickness of surface grid cell
        Δh = hⁿ[i, j, 1] - h⁻[i, j, 1]

        # Update surface salinity flux.
        # Note: the Δt below is the ocean time-step, eg.
        # ΔS = ⋯ - ∮ Qˢ dt ≈ ⋯ - Δtₒ * Qˢ 
        Qˢ[i, j, 1] += Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz])

        # Update previous ice thickness
        h⁻[i, j, 1] = hⁿ[i, j, 1]
    end
end

function ice_ocean_latent_heat!(coupled_model)
    ocean = coupled_model.ocean
    ice = coupled_model.ice
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    Qₒ = ice.model.external_thermal_fluxes.bottom
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    Δt = ocean.Δt
    hᵢ = ice.model.ice_thickness

    liquidus = ice.model.phase_transitions.liquidus
    grid = ocean.model.grid
    arch = architecture(grid)

    # What about the latent heat removed from the ocean when ice forms?
    # Is it immediately removed from the ocean? Or is it stored in the ice?
    launch!(arch, grid, :xy, _compute_ice_ocean_latent_heat!,
            Qₒ, grid, hᵢ, Tₒ, Sₒ, liquidus, ρₒ, cₒ, Δt)

    return nothing
end

@kernel function _compute_ice_ocean_latent_heat!(latent_heat,
                                                 grid,
                                                 ice_thickness,
                                                 ocean_temperature,
                                                 ocean_salinity,
                                                 liquidus,
                                                 ρₒ, cₒ, Δt)

    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)
    Qₒ = latent_heat
    hᵢ = ice_thickness
    Tₒ = ocean_temperature
    Sₒ = ocean_salinity

    δQ = zero(grid)
    icy_cell = @inbounds hᵢ[i, j, 1] > 0 # make ice bath approximation then

    @unroll for k = Nz:-1:1
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

        # Compute change in ocean thermal energy.
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
