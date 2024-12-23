using Oceananigans.Grids: inactive_node

""" Compute turbulent fluxes between an atmosphere and a surface state using similarity theory """
@kernel function _compute_atmosphere_surface_similarity_theory_fluxes!(fluxes_fields,
                                                                       coefficients,
                                                                       grid,
                                                                       clock,
                                                                       surface_state,
                                                                       surface_phase,
                                                                       surface_density,
                                                                       surface_heat_capacity,
                                                                       surface_temperature_units,
                                                                       surface_atmos_state,
                                                                       radiation,
                                                                       mole_fraction,
                                                                       vapor_saturation,
                                                                       atmosphere_reference_height,
                                                                       atmosphere_boundary_layer_height,
                                                                       atmos_thermodynamics_parameters)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) # index of the top ocean cell
    time = Time(clock.time)

    @inbounds begin
        uₐ = surface_atmos_state.u[i, j, 1]
        vₐ = surface_atmos_state.v[i, j, 1]
        Tₐ = surface_atmos_state.T[i, j, 1]
        pₐ = surface_atmos_state.p[i, j, 1]
        qₐ = surface_atmos_state.q[i, j, 1]
        Rs = surface_atmos_state.Qs[i, j, 1]
        Rℓ = surface_atmos_state.Qℓ[i, j, 1]

        # Extract state variables at cell centers
        # Ocean state
        uₒ = ℑxᶜᵃᵃ(i, j, kᴺ, grid, surface_state.u)
        vₒ = ℑyᵃᶜᵃ(i, j, kᴺ, grid, surface_state.v)
        Tₒ = surface_state.T[i, j, kᴺ]
        Tₒ = convert_to_kelvin(surface_temperature_units, Tₒ)
        Sₒ = surface_state.S[i, j, kᴺ]
    end

    # Build thermodynamic and dynamic states in the atmosphere and surface.
    # Notation:
    #   ⋅ 𝒬 ≡ thermodynamic state vector
    #   ⋅ 𝒰 ≡ "dynamic" state vector (thermodynamics + reference height + velocity)
    ℂₐ = atmos_thermodynamics_parameters
    𝒬ₐ = thermodynamic_atmospheric_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₐ, qₐ)

    hₐ = atmosphere_reference_height # elevation of atmos variables relative to surface
    
    Uₐ = SVector(uₐ, vₐ)
    𝒰ₐ = dynamic_atmos_state = SurfaceFluxes.StateValues(hₐ, Uₐ, 𝒬ₐ)

    # Build surface state with saturated specific humidity
    qₒ = seawater_saturation_specific_humidity(ℂₐ, Tₒ, Sₒ, 𝒬ₐ,
                                               mole_fraction,
                                               vapor_saturation,
                                               surface_phase)

    # Thermodynamic and dynamic (ocean) surface state:
    Uₒ = SVector(uₒ, vₒ)
     
    𝒬₀ = thermodynamic_surface_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₒ, qₒ)
    h₀ = zero(grid) # surface height
    𝒰₀ = dynamic_surface_state = SurfaceFluxes.StateValues(h₀, Uₒ, 𝒬₀)

    # Some parameters
    g = default_gravitational_acceleration
    
    surface_salinity = Sₒ
    prescribed_heat_fluxes = net_downwelling_radiation(i, j, grid, time, radiation, Rs, Rℓ) 
    radiative_properties = local_radiation_properties(i, j, kᴺ, grid, time, radiation)
    inactive_cell = inactive_node(i, j, kᴺ, grid, Center(), Center(), Center())

    turbulent_fluxes, surface_temperature = compute_turbulent_fluxes(coefficients,
                                                                     dynamic_surface_state, 
                                                                     dynamic_atmos_state, 
                                                                     prescribed_heat_fluxes,
                                                                     radiative_properties,
                                                                     surface_phase,
                                                                     surface_salinity,
                                                                     surface_density,
                                                                     surface_heat_capacity,
                                                                     mole_fraction,
                                                                     vapor_saturation,
                                                                     atmosphere_boundary_layer_height,
                                                                     ℂₐ, g, inactive_cell)

    # Store fluxes
    Qv = fluxes_fields.latent_heat
    Qc = fluxes_fields.sensible_heat
    Fv = fluxes_fields.water_vapor
    ρτx = fluxes_fields.x_momentum
    ρτy = fluxes_fields.y_momentum
    Ts  = fluxes_fields.T_surface

    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1]  = ifelse(inactive_cell, 0, turbulent_fluxes.latent_heat)
        Qc[i, j, 1]  = ifelse(inactive_cell, 0, turbulent_fluxes.sensible_heat)
        Fv[i, j, 1]  = ifelse(inactive_cell, 0, turbulent_fluxes.water_vapor)
        ρτx[i, j, 1] = ifelse(inactive_cell, 0, turbulent_fluxes.x_momentum)
        ρτy[i, j, 1] = ifelse(inactive_cell, 0, turbulent_fluxes.y_momentum)
        Ts[i, j, 1]  = surface_temperature
    end
end

@kernel function _assemble_atmosphere_sea_ice_fluxes!(centered_velocity_fluxes,
                                                      net_tracer_fluxes,
                                                      grid,
                                                      clock,
                                                      sea_ice_temperature,
                                                      sea_ice_salinity,
                                                      sea_ice_temperature_units,
                                                      turbulent_fluxes,
                                                      downwelling_radiation,
                                                      prescribed_freshwater_flux,
                                                      radiation_properties,
                                                      sea_ice_reference_density,
                                                      sea_ice_heat_capacity,
                                                      freshwater_density)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)

    @inbounds begin
        Tₒ = sea_ice_temperature[i, j, kᴺ]
        Tₒ = convert_to_kelvin(sea_ice_temperature_units, Tₒ)
        Sₒ = sea_ice_salinity[i, j, kᴺ]

        Qs = downwelling_radiation.shortwave[i, j, 1]
        Qℓ = downwelling_radiation.longwave[i, j, 1]

        Mp = prescribed_freshwater_flux[i, j, 1]

        Qc  = turbulent_fluxes.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv  = turbulent_fluxes.latent_heat[i, j, 1]   # latent heat flux
        Mv  = turbulent_fluxes.water_vapor[i, j, 1]   # mass flux of water vapor
        ρτx = turbulent_fluxes.x_momentum[i, j, 1]    # zonal momentum flux
        ρτy = turbulent_fluxes.y_momentum[i, j, 1]    # meridional momentum flux
    end

    # Compute heat fluxes, bulk flux first
    Qd = net_downwelling_radiation(i, j, grid, time, radiation_properties, Qs, Qℓ)
    Qu = net_upwelling_radiation(i, j, grid, time, radiation_properties, Tₒ)

    ΣQ = Qd + Qu + Qc + Qv

    # Convert from a mass flux to a volume flux (aka velocity)
    # by dividing with the density of freshwater.
    # Also switch the sign, for some reason we are given freshwater flux as positive down.
    ρf⁻¹ = 1 / freshwater_density
    ΣF   = - Mp * ρf⁻¹

    # Add the contribution from the turbulent water vapor flux, which has
    # a different sign convention as the prescribed water mass fluxes (positive upwards)
    Fv = Mv * ρf⁻¹
    ΣF += Fv

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    τx = centered_velocity_fluxes.u
    τy = centered_velocity_fluxes.v
    Jᵀ = net_tracer_fluxes.T
    Jˢ = net_tracer_fluxes.S

    ρₒ⁻¹ = 1 / sea_ice_reference_density
    cₒ   = sea_ice_heat_capacity

    atmos_ocean_τx = ρτx * ρₒ⁻¹
    atmos_ocean_τy = ρτy * ρₒ⁻¹
    atmos_ocean_Jᵀ = ΣQ  * ρₒ⁻¹ / cₒ
    atmos_ocean_Jˢ = - Sₒ * ΣF

    # Mask fluxes over land for convenience
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds begin
        τx[i, j, 1] = ifelse(inactive, 0, atmos_ocean_τx)
        τy[i, j, 1] = ifelse(inactive, 0, atmos_ocean_τy)
        Jᵀ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive, 0, atmos_ocean_Jˢ)
    end
end
