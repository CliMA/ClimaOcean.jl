""" The exchange fluxes depend on the atmosphere velocity but not the surface velocity """
struct WindVelocity end

""" The exchange fluxes depend on the relative velocity between the atmosphere and the surface """
struct RelativeVelocity end

"""
    buoyancy_scale(θ★, q★, 𝒬, ℂ, g)

Return the characteristic buoyancy scale `b★` associated with
the characteristic temperature `θ★`, specific humidity scale `q★`,
near-surface atmospheric thermodynamic state `𝒬`, thermodynamic
parameters `ℂ`, and gravitational acceleration `g`.

The buoyancy scale is defined in terms of the surface buoyancy flux,

```math
u★ b★ ≡ w′b′,
```

where `u★` is the friction velocity.
Using the definition of buoyancy for non-condensing air, we find that

```math
b★ = g / 𝒯ₐ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₐ * q★),
```
where ``𝒯ₐ`` is the virtual temperature of the atmosphere near the surface,
and ``δ = Rᵥ / R_d - 1``, where ``Rᵥ`` is the molar mass of water vapor and
``R_d`` is the molar mass of dry air.

Note that the Monin-Obukhov characteristic length scale is defined
in terms of `b★` and additionally the Von Karman constant `ϰ`,

```math
L★ = - u★² / ϰ b★ .
```
"""
@inline function buoyancy_scale(θ★, q★, 𝒬, ℂ, g)
    𝒯ₐ = AtmosphericThermodynamics.virtual_temperature(ℂ, 𝒬)
    qₐ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬)
    ε  = AtmosphericThermodynamics.Parameters.molmass_ratio(ℂ)
    δ  = ε - 1 # typically equal to 0.608

    b★ = g / 𝒯ₐ * (θ★ * (1 + δ * qₐ) + δ * 𝒯ₐ * q★)

    return b★
end

@inline velocity_differences(𝒰₁, 𝒰₀, ::RelativeVelocity) = @inbounds 𝒰₁.u[1] - 𝒰₀.u[1], 𝒰₁.u[2] - 𝒰₀.u[2]
@inline velocity_differences(𝒰₁, 𝒰₀, ::WindVelocity)     = @inbounds 𝒰₁.u[1], 𝒰₁.u[2] 

@inline function state_differences(ℂ, 𝒰₁, 𝒰₀, θ₀, S₀, Σ★, g, ρₒ, cpₒ, 
                                   water_mole_fraction,
                                   water_vapor_saturation,
                                   surface_temperature_type, 
                                   prescribed_heat_fluxes,
                                   radiative_properties,
                                   bulk_velocity,
                                   surface_phase)
    z₁ = 𝒰₁.z
    z₀ = 𝒰₀.z
    Δh = z₁ - z₀
    Δu, Δv = velocity_differences(𝒰₁, 𝒰₀, bulk_velocity)
    
    # Thermodynamic state
    𝒬₁ = 𝒰₁.ts
    𝒬₀ = 𝒰₀.ts

    ρₐ = AtmosphericThermodynamics.air_density(ℂ, 𝒬₁)
    cₚ = AtmosphericThermodynamics.cp_m(ℂ, 𝒬₁) # moist heat capacity
    ℰv = AtmosphericThermodynamics.latent_heat_vapor(ℂ, 𝒬₁)

    θ₀ = compute_surface_temperature(surface_temperature_type, θ₀, ℂ, 𝒬₀, ρₐ, cₚ, ℰv, Σ★, ρₒ, cpₒ, g,
                                     prescribed_heat_fluxes, 
                                     radiative_properties)

    θ₁ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬₁)

    # Temperature difference including the ``lapse rate'' `α = g / cₚ`
    Δθ = θ₁ - θ₀ + g / cₚ * Δh

    q₁ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬₁)

    # Recomputing the saturation specific humidity at the surface based on the new temperature
    q₀ = seawater_saturation_specific_humidity(ℂ, θ₀, S₀, 𝒬₁,
                                               water_mole_fraction,
                                               water_vapor_saturation,
                                               surface_phase)
    
    𝒬ₛ = AtmosphericThermodynamics.PhaseEquil_pTq(ℂ, 𝒬₀.p, θ₀, q₀)
    q₀ = AtmosphericThermodynamics.vapor_specific_humidity(ℂ, 𝒬ₛ)

    Δq = q₁ - q₀

    return Δh, Δu, Δv, Δθ, Δq, θ₀
end
