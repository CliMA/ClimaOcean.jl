# Boundary-layer diffusivity within hbl: cubic shape function G(σ), C¹ matching
# to interior K at σ = 1, and nonlocal transport coefficient γ.

@inline function shape_function(σ, G1, dG1)
    a₁ = σ - 2
    a₂ = 3 - 2σ
    a₃ = σ - 1
    return a₁ + a₂ * G1 + a₃ * dG1
end

# K(σ) = hbl · w(σ) · σ · (1 + σ · G(σ)).
@inline boundary_layer_diffusivity(σ, hbl, w, G1, dG1) =
    hbl * w * σ * (one(σ) + σ * shape_function(σ, G1, dG1))

# G(1) = Kint / (hbl·w);   dG/dσ|₁ = -dKdz / w + Cˢᵗ·Bo·Kint/u★⁴ (stable only).
# dG1 clamped ≤ 0 to prevent K growing upward at the BL base.
@inline function matching_coefficients(hbl, Kint, dKdz, w, Bo, u★, p)
    FT  = typeof(hbl)
    G1  = Kint / max(hbl * w, FT(1e-30))
    f₁  = ifelse(Bo ≥ zero(FT), p.Cˢᵗ * Bo / max(u★^4, FT(1e-30)), zero(FT))
    dG1 = - dKdz / max(w, FT(1e-30)) + f₁ * Kint
    return G1, min(dG1, zero(FT))
end

# Nonlocal-transport coefficient for tracers (zero for momentum and stable
# forcing). Derives MITgcm's `cg = C★·κᵥ·(Bˢ·κᵥ·ε)^(1/3) ≈ 6.33` and clamps
# γ at 100 s/m² to bound the runaway when ws·hbl is small.
@inline function nonlocal_transport(hbl, ws, Bo, p)
    FT = typeof(hbl)
    cg = p.C★ * p.κᵥ * cbrt(p.Bˢ * p.κᵥ * p.ε)
    γ  = min(cg / max(ws * hbl, FT(1e-30)), FT(100))
    return ifelse(Bo < zero(FT), γ, zero(FT))
end
