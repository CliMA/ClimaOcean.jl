# Interior viscosity / diffusivity at (Center, Center, Face) interfaces:
# IW background + shear-instability + convective-instability contributions.

@inline ϕ²(i, j, k, grid, ϕ, args...) = ϕ(i, j, k, grid, args...)^2

@inline function shear_squaredᶜᶜᶠ(i, j, k, grid, velocities)
    ∂z_u² = ℑxᶜᵃᵃ(i, j, k, grid, ϕ², ∂zᶠᶜᶠ, velocities.u)
    ∂z_v² = ℑyᵃᶜᵃ(i, j, k, grid, ϕ², ∂zᶜᶠᶠ, velocities.v)
    return ∂z_u² + ∂z_v²
end

# Smooth cubic shape: 1 at Ri ≤ 0, 0 at Ri ≥ Ri∞.
@inline function shear_factor(Ri, Ri∞, FT)
    r = min(max(Ri, zero(FT)) / Ri∞, one(FT))
    f = one(FT) - r * r
    return f * f * f
end

# Smooth cubic shape: 1 at N² ≤ N²ᶜᵒⁿ, 0 at N² ≥ 0. MITgcm uses a Heaviside
# switch instead; the smooth form is GPU-friendlier and otherwise matches.
@inline function convective_factor(N², N²ᶜᵒⁿ, FT)
    Ng = max(N², N²ᶜᵒⁿ)
    r  = min((N²ᶜᵒⁿ - Ng) / N²ᶜᵒⁿ, one(FT))
    f  = one(FT) - r * r
    return f * f * f
end

# IW background + shear + convective; matches MITgcm `kpp_routines.F:1197`.
# Returns (ν, κ) jointly so S², N², Ri and the smooth factors are computed once.
@inline function interior_diffusivitiesᶜᶜᶠ(i, j, k, grid, closure, velocities, tracers, buoyancy)
    FT = eltype(grid)
    p  = getclosure(i, j, closure).parameters

    S² = shear_squaredᶜᶜᶠ(i, j, k, grid, velocities)
    N² = ∂z_b(i, j, k, grid, buoyancy, tracers)
    Ri = N² / max(S², FT(1e-10))

    fSh = shear_factor(Ri, p.Ri∞, FT)
    fCv = convective_factor(N², p.N²ᶜᵒⁿ, FT)

    ν = p.νⁱʷ + fSh * p.ν₀ˢʰ + fCv * p.νᶜᵒⁿ
    κ = p.κⁱʷ + fSh * p.κ₀ˢʰ + fCv * p.κᶜᵒⁿ

    masked = peripheral_node(i, j, k, grid, Center(), Center(), Face())
    return ifelse(masked, zero(FT), ν), ifelse(masked, zero(FT), κ)
end
