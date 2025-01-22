@kernel function _assemble_atmosphere_sea_ice_fluxes!(net_heat_flux,
                                                      grid,
                                                      clock,
                                                      surface_temperature,
                                                      surface_temperature_units,
                                                      turbulent_fluxes,
                                                      downwelling_radiation,
                                                      stefan_boltzmann_constant,
                                                      albedo,
                                                      emissivity)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)

    @inbounds begin
        Tₒ = surface_temperature[i, j, kᴺ]
        Tₒ = convert_to_kelvin(surface_temperature_units, Tₒ)

        Qs  = downwelling_radiation.shortwave[i, j, 1]
        Qℓ  = downwelling_radiation.longwave[i, j, 1]
        Qc  = turbulent_fluxes.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv  = turbulent_fluxes.latent_heat[i, j, 1]   # latent heat flux
    end

    # Compute radiation fluxes
    σ = stefan_boltzmann_constant
    α = stateindex(albedo, i, j, 1, grid, time)
    ϵ = stateindex(emissivity, i, j, 1, grid, time)
    Qu = upwelling_radiation(Tₒ, σ, ϵ)
    Qd = net_downwelling_radiation(i, j, grid, time, α, ϵ, Qs, Qℓ)

    ΣQ = Qd + Qu + Qc + Qv

    # Mask fluxes over land for convenience
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds net_heat_flux[i, j, 1] = ifelse(inactive, zero(grid), ΣQ)
end

@kernel function reconstruct_momentum_fluxes!(grid, J, Jᶜᶜᶜ)
    i, j = @index(Global, NTuple)

    @inbounds begin
        J.u[i, j, 1] = ℑxᶠᵃᵃ(i, j, 1, grid, Jᶜᶜᶜ.u) 
        J.v[i, j, 1] = ℑyᵃᶠᵃ(i, j, 1, grid, Jᶜᶜᶜ.v) 
    end
end

@kernel function _assemble_atmosphere_ocean_fluxes!(net_fluxes,
                                                    grid,
                                                    clock,
                                                    interface_fluxes,
                                                    interface_temperature,
                                                    interior_state,
                                                    atmosphere_state,
                                                    interface_properties,
                                                    ocean_properties)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)
    time = Time(clock.time)

    @inbounds begin
        Tₛ = interface_temperature[i, j, kᴺ]
        Tₛ = convert_to_kelvin(ocean_properties.temperature_units, Tₛ)
        Sₒ = interior_state.S[i, j, kᴺ]

        Mp  = atmosphere_state.Mp[i, j, 1] # Prescribed freshwater flux
        Qs  = atmosphere_state.Qs[i, j, 1] # Downwelling shortwave radiation
        Qℓ  = atmosphere_state.Qℓ[i, j, 1] # Downwelling longwave radiation

        Qc  = interface_fluxes.sensible_heat[i, j, 1] # sensible or "conductive" heat flux
        Qv  = interface_fluxes.latent_heat[i, j, 1]   # latent heat flux
        Mv  = interface_fluxes.water_vapor[i, j, 1]   # mass flux of water vapor
        ρτx = interface_fluxes.x_momentum[i, j, 1]    # zonal momentum flux
        ρτy = interface_fluxes.y_momentum[i, j, 1]    # meridional momentum flux
    end

    radiation = interface_properties.radiation

    # Compute radiation fluxes
    σ = radiation.σ
    α = stateindex(radiation.α, i, j, 1, grid, time)
    ϵ = stateindex(radiation.ϵ, i, j, 1, grid, time)
    Qu = upwelling_radiation(Tₛ, σ, ϵ)
    Qd = net_downwelling_radiation(i, j, grid, time, α, ϵ, Qs, Qℓ)

    ΣQ = Qd + Qu + Qc + Qv

    # Convert from a mass flux to a volume flux (aka velocity)
    # by dividing with the density of freshwater.
    # Also switch the sign, for some reason we are given freshwater flux as positive down.
    ρf⁻¹ = 1 / ocean_properties.freshwater_density
    ΣF   = - Mp * ρf⁻¹

    # Add the contribution from the turbulent water vapor flux, which has
    # a different sign convention as the prescribed water mass fluxes (positive upwards)
    Fv = Mv * ρf⁻¹
    ΣF += Fv

    # Compute fluxes for u, v, T, S from momentum, heat, and freshwater fluxes
    τx = net_fluxes.u
    τy = net_fluxes.v
    Jᵀ = net_fluxes.T
    Jˢ = net_fluxes.S

    ρₒ⁻¹ = 1 / ocean_properties.reference_density
    cₒ   = ocean_properties.heat_capacity

    _τx = ρτx * ρₒ⁻¹
    _τy = ρτy * ρₒ⁻¹
    _Jᵀ = ΣQ  * ρₒ⁻¹ / cₒ
    _Jˢ = - Sₒ * ΣF

    # Mask fluxes over land for convenience
    inactive = inactive_node(i, j, kᴺ, grid, c, c, c)

    @inbounds begin
        τx[i, j, 1] = ifelse(inactive, zero(grid), _τx)
        τy[i, j, 1] = ifelse(inactive, zero(grid), _τy)
        Jᵀ[i, j, 1] = ifelse(inactive, zero(grid), _Jᵀ)
        Jˢ[i, j, 1] = ifelse(inactive, zero(grid), _Jˢ)
    end
end

