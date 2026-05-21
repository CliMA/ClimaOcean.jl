# Three-phase driver:
#   (1) cache interior (ν, κ) at every face,
#   (2) per-column scalars: hbl, u★, Bo, α, and matching coefs G1{u,s} / dG1{u,s},
#   (3) per-interface κu, κc, γ — pure lookup + per-face physics.

function compute_closure_fields!(diffusivities, closure::FlavorOfKPP, model; parameters = :xyz)
    arch  = model.architecture
    grid  = model.grid
    clock = model.clock

    radiation        = get_radiative_forcing(model)
    coriolis         = model.coriolis
    top_velocity_bcs = (u = model.velocities.u.boundary_conditions.top,
                        v = model.velocities.v.boundary_conditions.top)
    top_bcs = KPPTopBoundaryConditions(top_velocity_bcs, diffusivities.top_tracer_bcs.bcs)

    launch!(arch, grid, :xy, compute_kpp_column_fields!,
            diffusivities, grid, closure,
            model.velocities, model.tracers, model.buoyancy,
            top_bcs, radiation, coriolis, clock)

    launch!(arch, grid, parameters, compute_kpp_diffusivities!,
            diffusivities, grid, closure,
            model.buoyancy, radiation)

    return nothing
end

#####
##### Phase 2: column-level scalars (hbl, u★, Bo, α, matching coefficients)
#####

@kernel function compute_kpp_column_fields!(K, grid, closure, velocities, tracers, buoyancy,
                                            top_bcs, radiation, coriolis, clock)
    i, j = @index(Global, NTuple)

    FT = eltype(grid)
    Nz = grid.Nz
    p  = getclosure(i, j, closure).parameters
    fields = merge(velocities, tracers)

    u★  = friction_velocity(i, j, grid, clock, fields, top_bcs.velocities, p)
    Bo  = non_solar_buoyancy(i, j, grid, clock, fields, buoyancy, top_bcs.tracers)
    α   = αᶜᶜᶜ(i, j, grid, buoyancy, tracers)
    g   = buoyancy.formulation.gravitational_acceleration

    hbl = compute_boundary_layer_depth(i, j, grid, closure,
                                       velocities, tracers, buoyancy,
                                       u★, Bo, α, g, radiation, coriolis)

    # Capture cached (ν, κ) at the deepest face below hbl (subscript ₋) and the
    # first face above (subscript ₊) for the FD derivative dK/dz at hbl.
    z₀ = znode(i, j, Nz, grid, Center(), Center(), Center())
    ν₋ = zero(FT); ν₊ = zero(FT)
    κ₋ = zero(FT); κ₊ = zero(FT)
    z₋ = zero(FT); z₊ = zero(FT)
    crossed    = false
    have_below = false
    for k in 1:(Nz + 1)
        zf = znode(i, j, k, grid, Center(), Center(), Face())
        d  = z₀ - zf
        νₖ, κₖ = interior_diffusivitiesᶜᶜᶠ(i, j, k, grid, closure, velocities, tracers, buoyancy)
        @inbounds K.νᵢ[i, j, k] = νₖ
        @inbounds K.κᵢ[i, j, k] = κₖ

        below   = d > hbl
        capture = !below & !crossed

        ν₋ = ifelse(below,   νₖ, ν₋)
        κ₋ = ifelse(below,   κₖ, κ₋)
        z₋ = ifelse(below,   zf, z₋)
        ν₊ = ifelse(capture, νₖ, ν₊)
        κ₊ = ifelse(capture, κₖ, κ₊)
        z₊ = ifelse(capture, zf, z₊)

        crossed    = crossed    | !below
        have_below = have_below | below
    end

    # When hbl extends to the bottom, the FD pair is degenerate → dKdz = 0
    # so the matching reduces to a smooth G(σ) without spurious gradient terms.
    Δz  = max(z₊ - z₋, FT(1e-10))
    νh  = ν₋
    κh  = κ₋
    dνh = ifelse(have_below, (ν₊ - ν₋) / Δz, zero(FT))
    dκh = ifelse(have_below, (κ₊ - κ₋) / Δz, zero(FT))

    # Matching at σ = 1: column-level, so computed here once and reused per face.
    σ₁        = ifelse(Bo ≥ zero(FT), one(FT), p.ε)
    wm₁, ws₁  = velocity_scales(σ₁, hbl, u★, Bo, p)
    G1u, dG1u = matching_coefficients(hbl, νh, dνh, wm₁, Bo, u★, p)
    G1s, dG1s = matching_coefficients(hbl, κh, dκh, ws₁, Bo, u★, p)

    # Branchless land-column mask: zero everywhere on fully-land columns.
    wet = static_column_depthᶜᶜᵃ(i, j, grid) > zero(FT)
    @inbounds K.hbl[i, j, 1]  = ifelse(wet, hbl,  zero(FT))
    @inbounds K.u★[i, j, 1]   = ifelse(wet, u★,   zero(FT))
    @inbounds K.Bo[i, j, 1]   = ifelse(wet, Bo,   zero(FT))
    @inbounds K.α[i, j, 1]    = ifelse(wet, α,    zero(FT))
    @inbounds K.G1u[i, j, 1]  = ifelse(wet, G1u,  zero(FT))
    @inbounds K.dG1u[i, j, 1] = ifelse(wet, dG1u, zero(FT))
    @inbounds K.G1s[i, j, 1]  = ifelse(wet, G1s,  zero(FT))
    @inbounds K.dG1s[i, j, 1] = ifelse(wet, dG1s, zero(FT))
end

#####
##### Phase 3: per-interface κu, κc, γ
#####

@kernel function compute_kpp_diffusivities!(K, grid, closure, buoyancy, radiation)
    i, j, k = @index(Global, NTuple)
    _kpp_interface!(i, j, k, K, grid, closure, buoyancy, radiation)
end

@inline function _kpp_interface!(i, j, k, K, grid, closure, buoyancy, radiation)
    FT  = eltype(grid)
    Nz  = grid.Nz
    p   = getclosure(i, j, closure).parameters
    clo = getclosure(i, j, closure)

    @inbounds hbl  = K.hbl[i, j, 1]
    @inbounds u★   = K.u★[i, j, 1]
    @inbounds Bo   = K.Bo[i, j, 1]
    @inbounds α    = K.α[i, j, 1]
    @inbounds G1u  = K.G1u[i, j, 1]
    @inbounds dG1u = K.dG1u[i, j, 1]
    @inbounds G1s  = K.G1s[i, j, 1]
    @inbounds dG1s = K.dG1s[i, j, 1]

    @inbounds νᵢ = K.νᵢ[i, j, k]
    @inbounds κᵢ = K.κᵢ[i, j, k]

    g     = buoyancy.formulation.gravitational_acceleration
    z₀    = znode(i, j, Nz, grid, Center(), Center(), Center())
    d     = z₀ - znode(i, j, k, grid, Center(), Center(), Face())
    σ     = d / max(hbl, FT(1e-10))
    in_BL = (σ < one(FT)) & (σ ≥ zero(FT))

    # Local turbulent scales at this interface (SW-aware Bf).
    Bf     = buoyancy_forcing_above(i, j, d, Bo, radiation, α, g)
    σw     = ifelse(Bf ≥ zero(FT), one(FT), p.ε)
    wm, ws = velocity_scales(σw, hbl, u★, Bf, p)

    νᵇ = boundary_layer_diffusivity(σ, hbl, wm, G1u, dG1u)
    κᵇ = boundary_layer_diffusivity(σ, hbl, ws, G1s, dG1s)

    ν = min(ifelse(in_BL, max(νᵇ, νᵢ), νᵢ), clo.maximum_viscosity)
    κ = min(ifelse(in_BL, max(κᵇ, κᵢ), κᵢ), clo.maximum_diffusivity)
    γ = ifelse(in_BL, nonlocal_transport(hbl, ws, Bo, p), zero(FT))

    @inbounds K.κu[i, j, k] = ν
    @inbounds K.κc[i, j, k] = κ
    @inbounds K.γ[i, j, k]  = γ
end
