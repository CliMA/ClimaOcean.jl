# Boundary-layer depth `hbl` from the bulk-Richardson criterion (Large 1994
# eq. 21–23) plus Ekman / Monin-Obukhov clamps under stable forcing.

@inline function apply_stable_hbl_limits(hbl, u★, Bo, f, p)
    FT = typeof(hbl)
    hₑ = p.Cᵉᵏ * u★   / max(abs(f),    FT(1e-10))
    hₘ = p.Cᴹᴼ * u★^3 / max(p.κᵥ * Bo, FT(1e-10))
    return ifelse((Bo > zero(FT)) & p.limit_hbl_stable, min(hbl, hₑ, hₘ), hbl)
end

# Bulk-Ri sweep returning hbl. Branchless: sweep the whole column tracking
# the first crossing of Rib = Riᶜ via a `found::Bool` mask + ifelse.
@inline function compute_boundary_layer_depth(i, j, grid, closure,
                                              velocities, tracers, buoyancy,
                                              u★, Bo, α, g,
                                              radiation, coriolis)
    FT = eltype(grid)
    Nz = grid.Nz
    p  = getclosure(i, j, closure).parameters

    # Vt² coefficient (MITgcm `kpp_init_fixed.F:125`).
    βT  = FT(2//10)
    Vtc = p.Cᶜᵛ * sqrt(βT / (p.Bˢ * p.ε)) / (p.Riᶜ * p.κᵥ^2)

    # Surface reference at the top center.
    b₀ = buoyancy_perturbationᶜᶜᶜ(i, j, Nz, grid, buoyancy.formulation, tracers)
    u₀ = ℑxᶜᵃᵃ(i, j, Nz, grid, velocities.u)
    v₀ = ℑyᵃᶜᵃ(i, j, Nz, grid, velocities.v)
    z₀ = znode(i, j, Nz, grid, Center(), Center(), Center())
    H  = static_column_depthᶜᶜᵃ(i, j, grid)

    hbl   = H
    found = false
    Rib′  = zero(FT)
    d′    = zero(FT)

    for k in (Nz - 1):-1:1
        d   = z₀ - znode(i, j, k, grid, Center(), Center(), Center())
        Δb  = b₀ - buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, buoyancy.formulation, tracers)
        ΔV² = (u₀ - ℑxᶜᵃᵃ(i, j, k, grid, velocities.u))^2 +
              (v₀ - ℑyᵃᶜᵃ(i, j, k, grid, velocities.v))^2

        N²    = ℑzᵃᵃᶜ(i, j, k, grid, ∂z_b, buoyancy, tracers)
        Bf    = buoyancy_forcing_above(i, j, d, Bo, radiation, α, g)
        σ     = ifelse(Bf ≥ zero(FT), one(FT), p.ε)
        _, ws = velocity_scales(σ, d, u★, Bf, p)

        Vt² = d * ws * sqrt(max(N², zero(FT))) * Vtc
        Rib = d * Δb / max(ΔV² + Vt², FT(1e-10))
        Rib = ifelse(inactive_node(i, j, k, grid, Center(), Center(), Center()), zero(FT), Rib)

        crossed = (Rib ≥ p.Riᶜ) & !found
        hbl     = ifelse(crossed,
                         d′ + (d - d′) * (p.Riᶜ - Rib′) / max(Rib - Rib′, FT(1e-10)),
                         hbl)
        found   = found | (Rib ≥ p.Riᶜ)
        Rib′    = Rib
        d′      = d
    end

    f   = ℑxyᶜᶜᵃ(i, j, Nz, grid, fᶠᶠᵃ, coriolis)
    hbl = apply_stable_hbl_limits(hbl, u★, Bo, f, p)
    hbl = ifelse(found, hbl, p.minimum_boundary_layer_depth)
    return max(hbl, p.minimum_boundary_layer_depth)
end
