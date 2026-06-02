# Nonlocal-tracer-flux override. Tracer flux at face (i, j, k):
#
#     F = - κ · ∂c/∂z  -  κ · γ · Q₀
#         ╰──────────╯    ╰─────────╯
#         local            nonlocal (γ ≡ 0 for momentum and outside the BL)
#
# The nonlocal piece is independent of `c`, so it is treated explicitly even
# when the local piece runs through the vertically-implicit solver.

# Per-tracer Q₀ via integer-indexed lookup (Oceananigans passes `Val{tracer_index::Int}`).
@inline function surface_tracer_flux(i, j, grid, K, ::Val{id}, clock, fields) where id
    bc = K.top_tracer_bcs[id]
    return getbc(bc, i, j, grid, clock, fields)
end

@inline function nonlocal_z_flux(i, j, k, grid, K, val_id, clock, fields)
    FT = eltype(grid)
    κ  = @inbounds K.κc[i, j, k]
    γ  = @inbounds K.γ[i, j, k]
    Q₀ = surface_tracer_flux(i, j, grid, K, val_id, clock, fields)
    return ifelse(peripheral_node(i, j, k, grid, Center(), Center(), Face()),
                  zero(FT),
                  - κ * γ * Q₀)
end

# Explicit time discretization: full local + nonlocal flux.
@inline function diffusive_flux_z(i, j, k, grid, ::KPPVD, K, val_id::Val, c, clock, fields, buoyancy)
    κ = @inbounds K.κc[i, j, k]
    return - κ * ∂zᶜᶜᶠ(i, j, k, grid, c) + nonlocal_z_flux(i, j, k, grid, K, val_id, clock, fields)
end

# VITD on a periodic-z grid: local flux runs through the implicit solve;
# only the nonlocal piece is explicit.
@inline diffusive_flux_z(i, j, k, grid, ::VITD, ::KPPVD, K, val_id::Val, c, clock, fields, buoyancy) =
    nonlocal_z_flux(i, j, k, grid, K, val_id, clock, fields)

# VITD on a vertically-bounded grid: explicit (local + nonlocal) at the z
# boundaries, nonlocal-only in the bulk.
@inline function diffusive_flux_z(i, j, k, grid::VerticallyBoundedGrid, ::VITD, closure::KPPVD,
                                  K, val_id::Val, c, clock, fields, buoyancy)
    return ifelse((k == 1) | (k == grid.Nz + 1),
                  diffusive_flux_z(i, j, k, grid, ExplicitTimeDiscretization(), closure, K, val_id, c, clock, fields, buoyancy),
                  nonlocal_z_flux(i, j, k, grid, K, val_id, clock, fields))
end
