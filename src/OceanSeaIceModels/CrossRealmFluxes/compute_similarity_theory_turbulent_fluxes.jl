# See SurfaceFluxes.jl for other parameter set options.
default_businger_parameters(FT=Float64) = BusingerParams{FT}(Pr_0 = convert(FT, 0.74),
                                                             a_m  = convert(FT, 4.7),
                                                             a_h  = convert(FT, 4.7),
                                                             ζ_a  = convert(FT, 2.5),
                                                             γ    = convert(FT, 4.42))

@inline function compute_turbulent_surface_fluxes(similarity_function::BusingerParams{FT},
                                                  ::Nothing,
                                                  turbulent_fluxes,
                                                  atmos_state,
                                                  ocean_state) where FT

    # TODO: derive these from turbulent_fluxes.roughness_lengths
    zᵐ = zʰ = convert(FT, 5e-4) # τ = 0.3 => u★ = sqrt(τ / ρₐ) ~ z₀ ~ 5e-4

    # Solve for the surface fluxes with initial roughness length guess
    Uᵍ = zero(FT) # gustiness
    β = one(FT)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(atmos_state, ocean_state, zᵐ, zʰ, Uᵍ, β)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)

    fluxes = (;
        latent_heat         = conditions.lhf,
        sensible_heat       = conditions.shf,
        freshwater          = conditions.evaporation,
        zonal_momentum      = conditions.ρτxz,
        meridional_momentum = conditions.ρτyz
    )

    return fluxes
end

struct SimilarityFunction{FT}
    a :: FT
    b :: FT
    c :: FT
end

struct MomentumTracerSimilarityFunctions{U, C}
    momentum :: U
    tracers :: C
end


struct GravityWaveRoughnessLength{FT}
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
    air_kinematic_viscosity :: FT
end

struct AtmosphericState{Q, T, U, V}
    q :: Q
    θ :: T
    u :: U
    v :: V
end

AtmosphericState(q, θ, u) = AtmosphericState(q, θ, u, nothing)

@inline function (ψ::SimilarityFunction)(Ri)
    a = ψ.a
    b = ψ.b
    c = ψ.c

    ϕ⁻¹ = (1 - b * Ri)^c
    ψ_unstable = log((1 + ϕ⁻¹)^2 * (1 + ϕ⁻¹^2) / 8) - 2 * atan(ϕ⁻¹) + π/2
    ψ_stable = - a * Ri
    return ifelse(Ri < 0, ψ_unstable, ψ_stable)
end

@inline similarity_scale(ψ, h, ℓ, Ri) = 1 / (log(h/ℓ) - ψ(Ri) + ψ(ℓ * Ri / h))

function buoyancy_scale(θ★, q★, surface_state, parameters)
    θ★ = fluxes.θ
    q★ = fluxes.q
    𝒯₀ = virtual_temperature(parameters, surface_state)
    q₀ = surface_state.q
    θ₀ = surface_state.θ
    r = parameters.molar_mass_ratio
    g = parameters.gravitational_acceleration
    δ = r - 1
    b★ = g / 𝒯₀ * (θ★ * (1 + δ * q₀) + δ * θ₀ * q★)
    return b★
end

function fixed_point_fluxes(u★, θ★, q★,
                            surface_state,
                            inner_length_scales,
                            universal_function,
                            parameters)

    Δu = differences.u
    Δv = differences.v
    Δθ = differences.θ
    Δq = differences.q

    ϰ = parameters.von_karman_constant
    f = universal_function

    b★ = buoyancy_scale(θ★, q★, surface_state, parameters)
    Riₕ = - ϰ * h * b★ / u★^2

    ℓu = inner_length_scales.u(u★)
    ℓθ = inner_length_scales.θ(u★)
    ℓq = inner_length_scales.q(u★)

    χu = momentum_flux_scale(f, h, ℓu, Riₕ)
    χθ =   tracer_flux_scale(f, h, ℓθ, Riₕ)
    χq =   tracer_flux_scale(f, h, ℓq, Riₕ)

    u★ = ϰ * χu * sqrt(Δu^2 + Δv^2)
    θ★ = ϰ * χθ * Δθ
    q★ = ϰ * χq * Δq

    return u★, θ★, q★
end

function GravityWaveRoughnessLengths(FT=Float64;
                                     gravity_wave_parameter = 0.011,
                                     laminar_parameter = 0.11,
                                     air_kinematic_viscosity=1.5e-5)

    return GravityWaveRoughnessLengths(convert(FT, gravity_wave_parameter),
                                       convert(FT, laminar_parameter),
                                       convert(FT, air_kinematic_viscosity))
end

@inline function compute_turbulent_surface_fluxes(similarity_function::BusingerParams,
                                                  roughness_lengths::SimplifiedRoughnessLengths,
                                                  atmos_state,
                                                  ocean_state)

    # Solve for the surface fluxes with initial roughness length guess
    Uᵍ = zero(grid) # gustiness
    β = one(grid)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(atmos_state, ocean_State,
                                      roughness_lengths.momentum,
                                      roughness_lengths.heat
                                      Uᵍ, β)
    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)

    fluxes = (;
        latent_heat_flux         = conditions.lhf,
        sensible_heat_flux       = conditions.shf,
        freshwater_flux          = conditions.evaporation,
        zonal_momentum_flux      = conditions.ρτxz,
        meridional_momentum_flux = conditions.ρτyz,
    )

    return fluxes
end

@inline function compute_turbulent_surface_fluxes(roughness_lengths::GravityWaveRoughnessLengths,
                                                  atmos_state,
                                                  ocean_state)

    # Solve for the surface fluxes with initial roughness length guess
    Uᵍ = zero(grid) # gustiness
    β = one(grid)   # surface "resistance"
    values = SurfaceFluxes.ValuesOnly(atmos_state, ocean_State,
                                      roughness_lengths.momentum,
                                      roughness_lengths.heat
                                      Uᵍ, β)

    conditions = SurfaceFluxes.surface_conditions(turbulent_fluxes, values)

    fluxes = (;
        latent_heat_flux         = conditions.lhf,
        sensible_heat_flux       = conditions.shf,
        freshwater_flux          = conditions.evaporation,
        zonal_momentum_flux      = conditions.ρτxz,
        meridional_momentum_flux = conditions.ρτyz,
    )

    return fluxes
end

