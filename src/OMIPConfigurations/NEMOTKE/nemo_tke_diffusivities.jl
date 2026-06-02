# Diffusivity computation: Kₘ = K_ρ = Cᴷ · ℓ · √e, floored by background, capped by max.
# NEMO source: zdftke.F90 lines ~609-614.

@inline function nemo_eddy_coefficient(e_k, ℓ_k, p)
    FT = typeof(e_k)
    return p.Cᴷ * ℓ_k * sqrt(max(e_k, zero(FT)))
end

@inline viscosity_with_floors(K, p, K_max)   = min(max(K, p.νᵇ), K_max)
@inline diffusivity_with_floors(K, p, K_max) = min(max(K, p.κᵇ), K_max)
