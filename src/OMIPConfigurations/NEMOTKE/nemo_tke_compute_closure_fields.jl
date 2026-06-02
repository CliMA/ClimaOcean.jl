# NEMO TKE driver. Two Oceananigans hooks:
#
#   step_closure_prognostics!(K, closure, model, Δt)
#       — called by the timestepper with the substep Δt. Snapshot e → eⁿ at
#         substep 1, solve the column tridiagonal for eⁿ⁺¹, recompute ℓ.
#
#   compute_closure_fields!(K, closure, model)
#       — called by update_state! after the prognostic step. Computes the
#         per-face diffusivities Kₘ, K_ρ = Cᴷ · ℓ · √e with background
#         floors, viscosity/diffusivity caps, and EVD overwrite on N²<0.
#
# Tridiagonal scratch lives in K.γ, NOT K.eⁿ — eⁿ must persist across all
# RK3 substeps so each substep integrates one Δτ-step from the same eⁿ.
#
# Both kernels are written branchless: no `if/else` per cell, no `break` in
# loops, no early `return`. Per-cell masking uses `ifelse` so warps don't
# diverge across wet/dry columns or boundary/interior cells.

#####
##### Prognostic step: snapshot, TKE tridiagonal, mixing length.
#####

function step_closure_prognostics!(diffusivities, closure::FlavorOfNEMOTKE, model, Δt)
    arch  = model.architecture
    grid  = model.grid
    clock = model.clock

    if clock.stage == 1
        copyto!(parent(diffusivities.eⁿ), parent(diffusivities.e))
    end

    # Initial state construction passes Δt = Inf (no time step issued yet).
    isfinite(Δt) || return nothing

    velocities = model.velocities
    tracers    = model.tracers
    buoyancy   = model.buoyancy
    top_velocity_bcs = (u = model.velocities.u.boundary_conditions.top,
                        v = model.velocities.v.boundary_conditions.top)

    # Sea-ice fraction: 2D field if coupled, otherwise nothing → ℵ = 0 in the kernel.
    ℵ_field = nothing

    launch!(arch, grid, :xy, _compute_nemo_tke_column!,
            diffusivities, grid, closure, clock, Δt,
            velocities, tracers, buoyancy, top_velocity_bcs, ℵ_field)

    return nothing
end

#####
##### Per-face diffusivities only (called by update_state!).
#####

function compute_closure_fields!(diffusivities, closure::FlavorOfNEMOTKE, model; parameters = :xyz)
    arch     = model.architecture
    grid     = model.grid
    buoyancy = model.buoyancy
    tracers  = model.tracers

    launch!(arch, grid, parameters, _compute_nemo_tke_diffusivities!,
            diffusivities, grid, closure, buoyancy, tracers)

    return nothing
end

#####
##### Per-column kernel: cache surface forcing, solve TKE tridiagonal, build ℓ.
#####

@kernel function _compute_nemo_tke_column!(K, grid, closure, clock, Δt,
                                           velocities, tracers, buoyancy,
                                           top_velocity_bcs, ℵ_field)

    i, j = @index(Global, NTuple)
    FT   = eltype(grid)
    Nz   = grid.Nz
    p    = getclosure(i, j, closure).parameters

    fields = merge(velocities, tracers)

    # Surface-forcing cache (column scalars at k = 1).
    τx, τy = surface_stress_components(i, j, grid, clock, fields, top_velocity_bcs)
    u★²    = friction_velocity_squared(τx, τy)
    e_surf = surface_TKE(u★², p)
    ℵ      = _column_ice_fraction(i, j, ℵ_field, FT)

    @inbounds K.τx[i, j, 1]  = τx
    @inbounds K.τy[i, j, 1]  = τy
    @inbounds K.u★²[i, j, 1] = u★²
    @inbounds K.ℵ[i, j, 1]   = ℵ

    # Column-scope sources: Langmuir + wave-penetration setup.
    τ_mag = sqrt(τx*τx + τy*τy)
    u_s   = ifelse(p.apply_langmuir_circulation, stokes_velocity(τ_mag, p), zero(FT))
    u_s²  = u_s * u_s
    φ_deg = _column_latitude_degrees(i, j, Nz, grid)
    h_τ   = wave_decay_length(φ_deg, p)
    h_LC  = _diagnose_langmuir_depth(i, j, grid, buoyancy, tracers, u_s², Nz, FT)

    # TKE tridiagonal (Thomas, in-place on K.e; K.γ holds γ').
    # Row k = 1..Nz, with k = Nz the topmost cell:
    #     a_k · eⁿ⁺¹_{k-1} + b_k · eⁿ⁺¹_k + c_k · eⁿ⁺¹_{k+1} = d_k
    # Surface (k = Nz): Dirichlet eⁿ⁺¹ = e_surf  → row replaced with identity.
    # Bottom  (k = 1):  no flux → a_1 = 0 (set via Kᵇ = 0).
    @inbounds for k in 1:Nz
        dry    = peripheral_node(i, j, k, grid, Center(), Center(), Center())
        is_top = k == Nz
        is_bot = k == 1

        eⁿ_k    = K.eⁿ[i, j, k]
        ℓ_prev  = max(K.ℓ[i, j, k], p.minimum_mixing_length)
        ω_k     = p.Cᴰ * sqrt(max(eⁿ_k, zero(FT))) / ℓ_prev

        # Geometric coefficients. At the column ends, Kᵇ/Kᵃ are zeroed so a_k/c_k
        # vanish; we still index a valid neighbour to keep loads in-bounds.
        kᵇ      = ifelse(is_bot, k,         k - 1)
        kᵃ      = ifelse(is_top, k,         k + 1)
        Δz_c    = Δzᶜᶜᶜ(i, j, k,  grid)
        Δzᵇ     = ifelse(is_bot, one(FT),   Δzᶜᶜᶜ(i, j, kᵇ, grid))
        Δzᵃ     = ifelse(is_top, one(FT),   Δzᶜᶜᶜ(i, j, kᵃ, grid))
        Kᵇ      = ifelse(is_bot, zero(FT),  @inbounds(K.κu[i, j, k]))
        Kᵃ      = ifelse(is_top, zero(FT),  @inbounds(K.κu[i, j, k + 1]))

        a_k     = -Δt * Kᵇ / (Δz_c * Δzᵇ)
        c_k     = -Δt * Kᵃ / (Δz_c * Δzᵃ)
        Kᶜ      = FT(0.5) * (Kᵇ + Kᵃ)

        S²_k    = _shear_squared_centered(i, j, k, grid, velocities, FT)
        N²_k    = _N²_centered(i, j, k, grid, buoyancy, tracers, FT)
        z_c     = -znode(i, j, k, grid, Center(), Center(), Center())
        LC_k    = ifelse(p.apply_langmuir_circulation,
                         langmuir_source(FT(z_c), h_LC, u_s, p), zero(FT))
        WP_k    = ifelse(p.apply_wave_penetration,
                         wave_penetration_source(FT(z_c), e_surf, h_τ, ℵ, p), zero(FT))

        # Interior row.
        d_int   = eⁿ_k + Δt * (Kᶜ * S²_k - Kᶜ * N²_k + LC_k + WP_k)
        b_int   = one(FT) + Δt * ω_k - a_k - c_k

        # Surface-Dirichlet override at k = Nz.
        a_k     = ifelse(is_top, zero(FT), a_k)
        b_k     = ifelse(is_top, one(FT),  b_int)
        c_k     = ifelse(is_top, zero(FT), c_k)
        d_k     = ifelse(is_top, e_surf,   d_int)

        # Forward sweep. At k = 1, γ_prev/β_prev are zeroed (a_1 = 0 makes their
        # contribution drop out anyway, but we still load a valid index).
        γ_prev  = ifelse(is_bot, zero(FT), @inbounds(K.γ[i, j, kᵇ]))
        β_prev  = ifelse(is_bot, zero(FT), @inbounds(K.e[i, j, kᵇ]))
        denom   = b_k - a_k * γ_prev
        denom   = ifelse(abs(denom) < FT(1e-30), FT(1e-30), denom)
        γ_new   = c_k / denom
        β_new   = (d_k - a_k * β_prev) / denom

        # Dry cells: identity row → backward sweep yields the floor.
        @inbounds K.γ[i, j, k] = ifelse(dry, zero(FT),       γ_new)
        @inbounds K.e[i, j, k] = ifelse(dry, p.minimum_TKE,  β_new)
    end

    # Backward sweep: e_k = β_k - γ_k · e_{k+1}, floored. Top row already holds e_surf.
    @inbounds K.e[i, j, Nz] = max(K.e[i, j, Nz], p.minimum_TKE)
    @inbounds for k in (Nz - 1):-1:1
        dry   = peripheral_node(i, j, k, grid, Center(), Center(), Center())
        γ_k   = K.γ[i, j, k]
        e_new = max(K.e[i, j, k] - γ_k * K.e[i, j, k + 1], p.minimum_TKE)
        K.e[i, j, k] = ifelse(dry, p.minimum_TKE, e_new)
    end

    # Mixing length: per-cell natural ℓ₀, then two-pass gradient limiter.
    @inbounds for k in 1:Nz
        dry  = peripheral_node(i, j, k, grid, Center(), Center(), Center())
        N²_k = max(_N²_centered(i, j, k, grid, buoyancy, tracers, FT), zero(FT))
        ℓ_new = natural_length_scale(K.e[i, j, k], N²_k, p)
        K.ℓ[i, j, k] = ifelse(dry, p.minimum_mixing_length, ℓ_new)
    end

    # Top→bottom limit (k = Nz is the surface; walk Nz-1 down to 1).
    @inbounds for k in (Nz - 1):-1:1
        dry_k     = peripheral_node(i, j, k,     grid, Center(), Center(), Center())
        dry_above = peripheral_node(i, j, k + 1, grid, Center(), Center(), Center())
        skip      = dry_k | dry_above
        ℓ_lim     = min(K.ℓ[i, j, k + 1] + Δzᶜᶜᶜ(i, j, k + 1, grid), K.ℓ[i, j, k])
        K.ℓ[i, j, k] = ifelse(skip, K.ℓ[i, j, k], ℓ_lim)
    end
    # Bottom→top limit.
    @inbounds for k in 2:Nz
        dry_k     = peripheral_node(i, j, k,     grid, Center(), Center(), Center())
        dry_below = peripheral_node(i, j, k - 1, grid, Center(), Center(), Center())
        skip      = dry_k | dry_below
        ℓ_lim     = min(K.ℓ[i, j, k - 1] + Δzᶜᶜᶜ(i, j, k - 1, grid), K.ℓ[i, j, k])
        K.ℓ[i, j, k] = ifelse(skip, K.ℓ[i, j, k], ℓ_lim)
    end
end

#####
##### Per-face kernel: cache N², compute Kₘ, K_ρ at faces with EVD overwrite.
#####

@kernel function _compute_nemo_tke_diffusivities!(K, grid, closure, buoyancy, tracers)
    i, j, k = @index(Global, NTuple)
    _nemo_tke_face!(i, j, k, K, grid, closure, buoyancy, tracers)
end

@inline function _nemo_tke_face!(i, j, k, K, grid, closure, buoyancy, tracers)
    FT  = eltype(grid)
    Nz  = grid.Nz
    clo = getclosure(i, j, closure)
    p   = clo.parameters

    inactive = peripheral_node(i, j, k, grid, Center(), Center(), Face())

    # NaN guard against uninitialised T/S during the very first compute call.
    N²_raw  = ∂z_b(i, j, k, grid, buoyancy, tracers)
    N²_face = ifelse(inactive | !isfinite(N²_raw), zero(FT), N²_raw)
    @inbounds K.N²[i, j, k] = N²_face

    # Branchless face interpolation: clamp to interior centers so the average
    # collapses to the single-cell value at the top/bottom faces.
    k_lo   = clamp(k - 1, 1, Nz)
    k_hi   = clamp(k,     1, Nz)
    e_face = FT(0.5) * (@inbounds(K.e[i, j, k_lo]) + @inbounds(K.e[i, j, k_hi]))
    ℓ_face = FT(0.5) * (@inbounds(K.ℓ[i, j, k_lo]) + @inbounds(K.ℓ[i, j, k_hi]))

    Kᵀ  = nemo_eddy_coefficient(e_face, ℓ_face, p)
    Kₘ  = viscosity_with_floors(Kᵀ,   p, clo.maximum_viscosity)
    K_ρ = diffusivity_with_floors(Kᵀ, p, clo.maximum_diffusivity)

    Kₘ  = evd_overwrite_momentum(Kₘ, N²_face, p)
    K_ρ = evd_overwrite(K_ρ,         N²_face, p)

    Kₘ  = ifelse(inactive, zero(FT), Kₘ)
    K_ρ = ifelse(inactive, zero(FT), K_ρ)

    @inbounds K.κu[i, j, k] = Kₘ
    @inbounds K.κc[i, j, k] = K_ρ
    return nothing
end

#####
##### Local helpers
#####

# Sea-ice fraction at a column. Coupled-mode runs pass a 2D field; pure-ocean
# OMIP runs pass nothing → ℵ = 0.
@inline _column_ice_fraction(i, j, ::Nothing, ::Type{FT}) where {FT} = zero(FT)
@inline _column_ice_fraction(i, j, ℵ_field, ::Type{FT}) where {FT}   = @inbounds ℵ_field[i, j, 1]

# Latitude in degrees at a column. Defined on geographic grids; falls back to
# zero so the closure runs in single-column tests on RectilinearGrid.
@inline _column_latitude_degrees(i, j, k, grid) = _column_latitude_degrees(i, j, k, grid, eltype(grid))
@inline _column_latitude_degrees(i, j, k, grid::Oceananigans.Grids.LatitudeLongitudeGrid, ::Type{FT}) where {FT} =
    φnode(i, j, k, grid, Center(), Center(), Center())
@inline _column_latitude_degrees(i, j, k, grid::Oceananigans.Grids.OrthogonalSphericalShellGrid, ::Type{FT}) where {FT} =
    φnode(i, j, k, grid, Center(), Center(), Center())
@inline _column_latitude_degrees(i, j, k, grid::Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid, ::Type{FT}) where {FT} =
    _column_latitude_degrees(i, j, k, grid.underlying_grid, FT)
@inline _column_latitude_degrees(i, j, k, grid, ::Type{FT}) where {FT} = zero(FT)

# Centered shear: (∂_z u)² + (∂_z v)², averaged from faces to the cell centre.
# NaN guard for the very first compute call when fields are uninitialised.
@inline function _shear_squared_centered(i, j, k, grid, velocities, ::Type{FT}) where {FT}
    ∂z_u² = ℑxᶜᵃᵃ(i, j, k, grid, _ϕ², ∂zᶠᶜᶠ, velocities.u)
    ∂z_v² = ℑyᵃᶜᵃ(i, j, k, grid, _ϕ², ∂zᶜᶠᶠ, velocities.v)
    s²    = ∂z_u² + ∂z_v²
    return ifelse(isfinite(s²), s², zero(FT))
end

@inline _ϕ²(i, j, k, grid, ϕ, args...) = ϕ(i, j, k, grid, args...)^2

# Centered N²: average ∂z_b at k+1/2 and k-1/2 to the cell centre. Uses clamped
# face indices at column ends so the average collapses to a single ∂z_b call.
@inline function _N²_centered(i, j, k, grid, buoyancy, tracers, ::Type{FT}) where {FT}
    Nz   = grid.Nz
    k_lo = clamp(k,     2, Nz)
    k_hi = clamp(k + 1, 2, Nz)
    n²   = FT(0.5) * (∂z_b(i, j, k_lo, grid, buoyancy, tracers) +
                      ∂z_b(i, j, k_hi, grid, buoyancy, tracers))
    return ifelse(isfinite(n²), n², zero(FT))
end

# Diagnostic Langmuir depth: deepest level where the cumulative buoyancy
# integral max(N², 0) · |z| · Δz exceeds u_s²/2. NEMO zdftke L341-369.
# Branchless: walk the entire column without `break`; `h_LC` records only the
# first crossing.
@inline function _diagnose_langmuir_depth(i, j, grid, buoyancy, tracers, u_s², Nz, ::Type{FT}) where {FT}
    threshold = u_s² * FT(0.5)
    cum       = zero(FT)
    h_LC      = zero(FT)
    hit       = false
    @inbounds for k in Nz:-1:1
        dry    = peripheral_node(i, j, k, grid, Center(), Center(), Center())
        N²_raw = ∂z_b(i, j, k, grid, buoyancy, tracers)
        N²_k   = ifelse(isfinite(N²_raw), max(N²_raw, zero(FT)), zero(FT))
        z_c    = -znode(i, j, k, grid, Center(), Center(), Center())
        Δz_k   = Δzᶜᶜᶜ(i, j, k, grid)
        cum   += ifelse(dry, zero(FT), N²_k * z_c * Δz_k)
        crossed   = cum > threshold
        first_hit = crossed & !hit
        h_LC      = ifelse(first_hit, z_c, h_LC)
        hit       = hit | crossed
    end
    return h_LC
end
