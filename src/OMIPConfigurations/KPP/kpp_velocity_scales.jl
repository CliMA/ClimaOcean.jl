# Monin-Obukhov velocity scales wm (momentum) and ws (scalars). Stable: linear
# in ζ. Unstable: positive-root form `κᵥ·u★·(...)^(±)`, robust at u★ → 0.

@inline function velocity_scales(σ, hbl, u★, Bf, params)
    FT  = typeof(σ)
    κᵥ  = params.κᵥ
    Cˢᵗ = params.Cˢᵗ

    ζ  = κᵥ * σ * hbl * Bf / max(u★^3, FT(1e-20))
    ζu = min(ζ, zero(FT))

    # Stable shared scale.
    w = κᵥ * u★ / (one(FT) + Cˢᵗ * max(ζ, zero(FT)))

    wm = ifelse(Bf ≥ zero(FT), w,
                κᵥ * u★ * ifelse(ζu > params.ζᵐ,
                                 sqrt(sqrt(one(FT) - params.Cᵐ * ζu)),
                                 cbrt(params.Aᵐ - params.Bᵐ * min(ζu, params.ζᵐ))))

    ws = ifelse(Bf ≥ zero(FT), w,
                κᵥ * u★ * ifelse(ζu > params.ζˢ,
                                 sqrt(one(FT) - params.Cˢ * ζu),
                                 cbrt(params.Aˢ - params.Bˢ * min(ζu, params.ζˢ))))

    return wm, ws
end
