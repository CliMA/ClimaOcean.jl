struct ConstantCoefficientFluxes{FT, U}
    momentum      :: FT
    heat          :: FT
    water_vapor   :: FT
    bulk_velocity :: U
end

@inline function compute_turbulent_fluxes(coefficients::ConstantCoefficientFluxes,
                                          surface_state,
                                          atmos_state,
                                          prescribed_heat_fluxes, # Possibly use in state_differences
                                          radiative_properties,
                                          surface_phase,
                                          surface_salinity,
                                          surface_density,
                                          surface_heat_capacity,
                                          mole_fraction,
                                          vapor_saturation,
                                          atmos_boundary_layer_height,
                                          thermodynamics_parameters,
                                          gravitational_acceleration,
                                          inactive_cell)
    
    # Prescribed difference between two states
    ℂₐ = thermodynamics_parameters
    FT = eltype(ℂₐ)

    # state differences...
    Δh, Δu, Δv, Δθ, Δq, θs = state_differences(thermodynamics_parameters,
                                               atmos_state,
                                               surface_state,
                                               surface_temperature,
                                               surface_salinity,
                                               nothing, # Σ★, required only for a `SkinTemperature` type (not needed here)
                                               gravitational_acceleration,
                                               surface_density,
                                               surface_heat_capacity,
                                               mole_fraction,
                                               vapor_saturation,
                                               BulkTemperature(), # There is no solver for this!
                                               prescribed_heat_fluxes,
                                               radiative_properties,
                                               coefficients.bulk_velocity,
                                               surface_phase)

    # The inital velocity scale assumes that the gustiness velocity `Uᴳ` is equal to 0.5 ms⁻¹. 
    # The initial surface temperature is the same as the surface temperature.
    # These will be refined later on.
    ΔU = sqrt(Δu^2 + Δv^2)
    
    # `u★² ≡ sqrt(τx² + τy²)`
    # We remove the gustiness by dividing by `ΔU`
    τx = - coefficients.momentum * Δu * ΔU
    τy = - coefficients.momentum * Δv * ΔU

    𝒬ₐ = atmos_state.ts
    ρₐ = AtmosphericThermodynamics.air_density(ℂₐ, 𝒬ₐ)
    cₚ = AtmosphericThermodynamics.cp_m(ℂₐ, 𝒬ₐ) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂₐ, 𝒬ₐ)

    fluxes = (;
        sensible_heat = - ρₐ * cₚ * ΔU * Δθ,
        latent_heat   = - ρₐ * ΔU * Δq * ℰv,
        water_vapor   = - ρₐ * ΔU * Δq,
        x_momentum    = + ρₐ * τx,
        y_momentum    = + ρₐ * τy,
    )
    
    return fluxes, θs
end
