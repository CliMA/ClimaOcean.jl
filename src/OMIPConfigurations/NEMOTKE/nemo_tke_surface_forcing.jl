# Surface forcing: τ extraction, u★² computation, Dirichlet surface BC for e.
# NEMO source: zdftke.F90 lines ~310-314 (surface BC).

# τ at a column: pull from the top boundary conditions of u and v. Mirrors
# the pattern used by KPP's friction_velocity (kpp_surface_forcing.jl).
@inline function surface_stress_components(i, j, grid, clock, fields, top_velocity_bcs)
    # τ in the same units NEMO uses: kinematic (m²/s²) — the top BCs of u/v in
    # Oceananigans return tendency in m²/s² already.
    τx = Oceananigans.BoundaryConditions.getbc(top_velocity_bcs.u, i, j, grid, clock, fields)
    τy = Oceananigans.BoundaryConditions.getbc(top_velocity_bcs.v, i, j, grid, clock, fields)
    return τx, τy
end

# u★² = |τ| (kinematic stress, m²/s²)
@inline function friction_velocity_squared(τx, τy)
    return sqrt(τx^2 + τy^2)
end

# NEMO surface TKE BC: e(top) = max(rn_emin0, rn_ebb · u★²)
@inline function surface_TKE(u★², p)
    return max(p.minimum_surface_TKE, p.Cᵇ * u★²)
end
