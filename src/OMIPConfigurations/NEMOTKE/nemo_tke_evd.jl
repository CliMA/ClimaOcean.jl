# Enhanced Vertical Diffusion (EVD), zdfevd.F90 with `nn_evdm=1` (apply to both).
# When N²(face) ≤ -1e-12, overwrite Kₘ and K_ρ to κᶜⁿᵛ (rn_avevd = 100 m²/s).
# Documented deviation from NEMO: NEMO uses min(rn2, rn2b); we use only N²
# (Oceananigans is single-step, no leapfrog history).

@inline function evd_overwrite(K_face, N²_face, p)
    FT = typeof(K_face)
    triggered = (N²_face <= FT(-1e-12)) & p.apply_enhanced_vertical_diffusion
    return ifelse(triggered, p.κᶜⁿᵛ, K_face)
end

@inline function evd_overwrite_momentum(K_face, N²_face, p)
    FT = typeof(K_face)
    triggered = (N²_face <= FT(-1e-12)) & p.apply_enhanced_vertical_diffusion & p.apply_evd_to_momentum
    return ifelse(triggered, p.κᶜⁿᵛ, K_face)
end
