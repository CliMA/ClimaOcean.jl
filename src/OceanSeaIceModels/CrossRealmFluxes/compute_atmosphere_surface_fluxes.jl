""" Compute turbulent fluxes between an atmosphere and a surface state using similarity theory """
@kernel function _compute_atmosphere_surface_similarity_theory_fluxes!(similarity_theory,
                                                                       grid,
                                                                       clock,
                                                                       surface_state,
                                                                       surface_density,
                                                                       surface_heat_capacity,
                                                                       surface_temperature_units,
                                                                       surface_atmos_state,
                                                                       radiation,
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
    surface_type = AtmosphericThermodynamics.Liquid()
    qₒ = seawater_saturation_specific_humidity(ℂₐ, Tₒ, Sₒ, 𝒬ₐ,
                                               similarity_theory.water_mole_fraction,
                                               similarity_theory.water_vapor_saturation,
                                               surface_type)

    # Thermodynamic and dynamic (ocean) surface state:
    Uₒ = SVector(uₒ, vₒ)
     
    𝒬₀ = thermodynamic_surface_state = AtmosphericThermodynamics.PhaseEquil_pTq(ℂₐ, pₐ, Tₒ, qₒ)
    h₀ = zero(grid) # surface height
    𝒰₀ = dynamic_surface_state = SurfaceFluxes.StateValues(h₀, Uₒ, 𝒬₀)

    # Some parameters
    g = default_gravitational_acceleration
    ϰ = similarity_theory.von_karman_constant
    
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)
    maxiter  = ifelse(inactive, 1, similarity_theory.maxiter)

    prescribed_heat_fluxes = net_downwelling_radiation(i, j, grid, time, radiation, Rs, Rℓ) 
    radiative_properties = local_radiation_properties(i, j, kᴺ, grid, time, radiation)

    turbulent_fluxes, surface_temperature = compute_similarity_theory_fluxes(similarity_theory,
                                                                             dynamic_surface_state, 
                                                                             dynamic_atmos_state, 
                                                                             prescribed_heat_fluxes,
                                                                             radiative_properties,
                                                                             Sₒ,
                                                                             surface_density,
                                                                             surface_heat_capacity,
                                                                             atmosphere_boundary_layer_height,
                                                                             ℂₐ, g, ϰ, maxiter)

    # Store fluxes
    Qv = similarity_theory.fields.latent_heat
    Qc = similarity_theory.fields.sensible_heat
    Fv = similarity_theory.fields.water_vapor
    ρτx = similarity_theory.fields.x_momentum
    ρτy = similarity_theory.fields.y_momentum
    Ts  = similarity_theory.fields.T_surface

    @inbounds begin
        # +0: cooling, -0: heating
        Qv[i, j, 1]  = ifelse(inactive, 0, turbulent_fluxes.latent_heat)
        Qc[i, j, 1]  = ifelse(inactive, 0, turbulent_fluxes.sensible_heat)
        Fv[i, j, 1]  = ifelse(inactive, 0, turbulent_fluxes.water_vapor)
        ρτx[i, j, 1] = ifelse(inactive, 0, turbulent_fluxes.x_momentum)
        ρτy[i, j, 1] = ifelse(inactive, 0, turbulent_fluxes.y_momentum)
        Ts[i, j, 1]  = surface_temperature
    end
end
