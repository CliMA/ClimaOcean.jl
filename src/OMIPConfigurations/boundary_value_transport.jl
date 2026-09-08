using KernelAbstractions: @index, @kernel
using Oceananigans.Architectures: architecture
using Oceananigans.BuoyancyFormulations: ∂z_b
using Oceananigans.Fields: VelocityFields
using Oceananigans.Grids: Center, Face, peripheral_node
using Oceananigans.Operators: Az⁻¹ᶜᶜᶠ, Δx_qᶜᶠᶠ, Δy_qᶠᶜᶠ, Δz⁻¹ᶜᶠᶜ, Δz⁻¹ᶠᶜᶜ,
                              Δzᶜᶠᶜ, Δzᶜᶠᶠ, Δzᶠᶜᶜ, Δzᶠᶜᶠ,
                              δxᶜᵃᵃ, δyᵃᶜᵃ, δzᵃᵃᶜ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Oceananigans.TurbulenceClosures: AbstractTurbulenceClosure, ExplicitTimeDiscretization,
                                       FluxTapering, convert_diffusivity, issd_coefficient_loc,
                                       κᶜᶠᶠ, κᶠᶜᶠ, ϵSxᶠᶜᶠ, ϵSyᶜᶠᶠ
using Oceananigans.Utils: KernelParameters, launch!, prettysummary
using Adapt: Adapt, adapt

# The mesoscale eddy transport of Ferrari, Griffies, Nurser & Vallis (2010), "A boundary-value
# problem for the parameterized mesoscale eddy transport", Ocean Modelling 32, 143--156.
#
# Instead of setting the eddy-induced transport to the local Gent et al. (1995) value Υᴳᴹ = -κ S,
# the transport in each column solves the two-point boundary-value problem (their Eq. 16)
#
#     (c² ∂²/∂z² - N²) Υ = (g/ρₒ) κ ∇_z ρ = -N² Υᴳᴹ ,      Υ(η) = Υ(-H) = 0 ,
#
# In modal space this is a low-pass filter, Υₘ = Υᴳᴹₘ / (1 + (c/cₘ)²), so the transport is dominated
# by the gravest baroclinic modes — the vertical structure geostrophic turbulence actually selects
# (Smith and Vallis 2002). Three consequences matter for a global run:
#
#   * the homogeneous Dirichlet conditions are satisfied by construction, so no tapering of κ and no
#     boundary-layer matching is needed at the surface or the bottom, where Υᴳᴹ = -κ S does not vanish;
#   * the second-order operator interpolates through weakly stratified layers, so neither a floor on
#     N² nor a ceiling on the neutral slope is required to regularize mixed layers or mode water;
#   * for the nonlinear Eady spindown the solution matches the Fox-Kemper et al. (2008) structure
#     function almost exactly (their Fig. 1), where Υᴳᴹ is a step function with delta-function jets.
#
# The transport enters the model as an eddy-induced velocity (u★, w★) = (∂z Υ, -∇ ⋅ Υ), the same
# advective route Oceananigans' `AdvectiveFormulation` takes; Redi mixing is untouched and stays in
# a companion `IsopycnalSkewSymmetricDiffusivity` carrying `κ_symmetric` alone.
#
# TODO: this belongs in Oceananigans as a third `skew_flux_formulation` alongside
# `DiffusiveFormulation` and `AdvectiveFormulation`, since it reuses that machinery wholesale
# (it only inserts a column solve between the transport and the eddy velocities). It lives here
# while we evaluate whether it fixes the AMOC.

"""
    struct BoundaryValueTransport{K, L, FT}

Parameterized mesoscale eddy transport determined by the boundary-value problem of Ferrari et al.
(2010), applied as an eddy-induced advection. Carries the skew diffusivity `κ_skew` only: pair it
with an `IsopycnalSkewSymmetricDiffusivity(κ_skew=nothing, κ_symmetric=...)` for Redi mixing.

`κ_skew` must be depth-independent — the boundary-value problem, not the diffusivity, supplies the
vertical structure of the transport (Ferrari et al. 2010, Section 4.1).
"""
struct BoundaryValueTransport{K, L, FT} <: AbstractTurbulenceClosure{ExplicitTimeDiscretization, 2}
    κ_skew        :: K
    slope_limiter :: L
    mode_number   :: FT   # M in c_M = (M π)⁻¹ ∫ N dz
    minimum_speed :: FT   # c_min
    minimum_N²    :: FT   # N₀², bounding N² from below so the scheme stays a sink of potential energy
end

"""
    BoundaryValueTransport(FT = Oceananigans.defaults.FloatType;
                           κ_skew,
                           slope_limiter = FluxTapering(1e-2),
                           mode_number = 2,
                           minimum_speed = 0.1,
                           minimum_N² = 1e-24)

Build the Ferrari et al. (2010) eddy transport closure.

The speed weighting the second-order operator is `c = max(minimum_speed, c_M)` with
`c_M = (M π)⁻¹ ∫ N dz` the WKB baroclinic gravity wave speed of mode `M` (their Eqs. 34 and 60).
`M = 1` gives the first baroclinic mode, whose amplitude is about half the truncated Gent et al.
(1995) transport; `M > 1` filters less and restores the amplitude while mode one still dominates
the vertical structure (their Section 4.5.1). The `minimum_speed` floor keeps the transport bounded
in weakly stratified columns and relaxes the `Δz < c/N` requirement the operator otherwise imposes.
"""
function BoundaryValueTransport(FT = Oceananigans.defaults.FloatType;
                                κ_skew,
                                slope_limiter = FluxTapering(1e-2),
                                mode_number = 2,
                                minimum_speed = 0.1,
                                minimum_N² = 1e-24)

    mode_number > 0 || throw(ArgumentError("mode_number must be positive, got $mode_number"))
    minimum_speed > 0 || throw(ArgumentError("minimum_speed must be positive, got $minimum_speed"))

    return BoundaryValueTransport(convert_diffusivity(FT, κ_skew),
                                  slope_limiter,
                                  convert(FT, mode_number),
                                  convert(FT, minimum_speed),
                                  convert(FT, minimum_N²))
end

Adapt.adapt_structure(to, closure::BoundaryValueTransport) =
    BoundaryValueTransport(adapt(to, closure.κ_skew),
                           adapt(to, closure.slope_limiter),
                           closure.mode_number,
                           closure.minimum_speed,
                           closure.minimum_N²)

Base.summary(::BoundaryValueTransport) = "BoundaryValueTransport"

function Base.show(io::IO, closure::BoundaryValueTransport)
    print(io, "BoundaryValueTransport (Ferrari et al. 2010)", '\n',
              "├── κ_skew: ",        prettysummary(closure.κ_skew), '\n',
              "├── slope_limiter: ", prettysummary(closure.slope_limiter), '\n',
              "├── mode_number: ",   prettysummary(closure.mode_number), '\n',
              "├── minimum_speed: ", prettysummary(closure.minimum_speed), '\n',
              "└── minimum_N²: ",    prettysummary(closure.minimum_N²))
end

#####
##### Closure fields
#####
#
# Storage layout:
#   u, v, w  — eddy-induced velocity, at the velocity locations, read by the tracer advection
#              through `closure_auxiliary_velocity`.
#   Υˣ, Υʸ   — the horizontal eddy transport at (Face, Center, Face) and (Center, Face, Face).
#              Holds Υᴳᴹ = -κ S after the first kernel and is overwritten in place with the
#              boundary-value solution by the second, so no extra column of storage is needed.
#   t        — scratch for the modified upper diagonal of the tridiagonal sweep. Only ever read at
#              the column it was written in, so its location is immaterial and both components
#              reuse it.

function Oceananigans.TurbulenceClosures.build_closure_fields(grid, clock, tracer_names, bcs,
                                                              closure::BoundaryValueTransport)
    velocities = VelocityFields(grid)

    return merge(velocities, (Υˣ = Field{Face, Center, Face}(grid),
                              Υʸ = Field{Center, Face, Face}(grid),
                              t  = Field{Center, Center, Face}(grid)))
end

@inline Oceananigans.TurbulenceClosures.closure_auxiliary_velocity(::BoundaryValueTransport, K, val_tracer_name) =
    (u = K.u, v = K.v, w = K.w)

function Oceananigans.TurbulenceClosures.compute_closure_fields!(closure_fields,
                                                                 closure::BoundaryValueTransport,
                                                                 model; parameters = :xyz)

    grid     = model.grid
    arch     = architecture(grid)
    clock    = model.clock
    buoyancy = Oceananigans.TurbulenceClosures.buoyancy_force(model)
    tracers  = Oceananigans.TurbulenceClosures.buoyancy_tracers(model)

    Nx, Ny, Nz = size(grid)

    # The eddy velocities take a horizontal difference of Υ, so Υ is needed one point into the halo.
    # Solving there too — every column is independent — avoids a halo exchange, and with it the
    # vector-component sign flip the tripolar fold would otherwise demand of a transport component.
    transport_parameters = KernelParameters((Nx + 2, Ny + 2, Nz + 1), (-1, -1, 0))
    column_parameters    = KernelParameters((Nx + 2, Ny + 2),         (-1, -1))

    launch!(arch, grid, transport_parameters, _compute_eddy_transport!,
            closure_fields.Υˣ, closure_fields.Υʸ, grid, clock, closure.κ_skew,
            closure.slope_limiter, buoyancy, fields(model))

    launch!(arch, grid, column_parameters, _solve_transport_boundary_value_problem!,
            closure_fields.Υˣ, closure_fields.Υʸ, closure_fields.t, grid, buoyancy, tracers,
            closure.mode_number, closure.minimum_speed, closure.minimum_N²)

    launch!(arch, grid, parameters, _compute_bvp_eddy_velocities!,
            closure_fields.u, closure_fields.v, closure_fields.w, grid,
            closure_fields.Υˣ, closure_fields.Υʸ)

    return nothing
end

#####
##### The local Gent et al. (1995) transport, Υᴳᴹ = -κ S
#####

@kernel function _compute_eddy_transport!(Υˣ, Υʸ, grid, clock, κ, slope_limiter, buoyancy, fields)
    i, j, k = @index(Global, NTuple)

    t  = clock.time
    κˣ = κᶠᶜᶠ(i, j, k, grid, issd_coefficient_loc, κ, t, fields)
    κʸ = κᶜᶠᶠ(i, j, k, grid, issd_coefficient_loc, κ, t, fields)

    Sˣ = ϵSxᶠᶜᶠ(i, j, k, grid, slope_limiter, buoyancy, fields)
    Sʸ = ϵSyᶜᶠᶠ(i, j, k, grid, slope_limiter, buoyancy, fields)

    # The slope stencil reaches into the halo, where the folded metrics of a tripolar grid can
    # return a non-finite slope; one such level would poison the whole column solve below.
    @inbounds begin
        Υˣ[i, j, k] = - finite_or_zero(κˣ * Sˣ, grid)
        Υʸ[i, j, k] = - finite_or_zero(κʸ * Sʸ, grid)
    end
end

#####
##### The boundary-value problem, solved column by column
#####

# N² at the two transport locations, interpolated from the (Center, Center, Face) buoyancy gradient.
@inline N²ᶠᶜᶠ(i, j, k, grid, buoyancy, tracers) = ℑxᶠᵃᵃ(i, j, k, grid, ∂z_b, buoyancy, tracers)
@inline N²ᶜᶠᶠ(i, j, k, grid, buoyancy, tracers) = ℑyᵃᶠᵃ(i, j, k, grid, ∂z_b, buoyancy, tracers)

# Squared baroclinic gravity wave speed weighting the second-order operator, c = max(c_min, c_M)
# with the WKB estimate c_M = (M π)⁻¹ ∫ N dz (Ferrari et al. 2010, Eqs. 34 and 60). N² is known at
# the interior faces — the unknowns of the problem — while the integral runs over the whole column,
# so the estimate is written as the column mean of N times the column depth. Summing N Δz over the
# interior faces alone would truncate the integral by one cell, which biases c at first order in Δz.
@inline function squared_wave_speed(i, j, grid, buoyancy, tracers, ℓx, ℓy, N², Δzᶜ, Δzᶠ, M, c★, N²₀)
    Nz = size(grid, 3)

    ∫Ndz = zero(grid)
    ∫dz  = zero(grid)
    depth = zero(grid)

    for k in 2:Nz
        wet = !peripheral_node(i, j, k, grid, ℓx, ℓy, Face())
        N²ᵏ = finite_or_zero(N²(i, j, k, grid, buoyancy, tracers), grid)
        Δz  = Δzᶠ(i, j, k, grid)
        ∫Ndz += ifelse(wet, sqrt(max(N²ᵏ, N²₀)) * Δz, zero(grid))
        ∫dz  += ifelse(wet, Δz, zero(grid))
    end

    for k in 1:Nz
        wet = !peripheral_node(i, j, k, grid, ℓx, ℓy, Center())
        depth += ifelse(wet, Δzᶜ(i, j, k, grid), zero(grid))
    end

    N̄ = ifelse(∫dz > 0, ∫Ndz / ∫dz, zero(grid))
    c = max(c★, N̄ * depth / (M * π))

    return c^2
end

# Thomas algorithm on one column. Row k of the discrete problem at an interior face reads
#
#     aₖ Υₖ₋₁ + bₖ Υₖ + cₖ Υₖ₊₁ = dₖ ,     dₖ = -N²ₖ Υᴳᴹₖ ,
#
# with aₖ = c²/(Δzᶜₖ₋₁ Δzᶠₖ), cₖ = c²/(Δzᶜₖ Δzᶠₖ) and bₖ = -(aₖ + cₖ) - N²ₖ. Faces on the boundary
# of the fluid — the surface, the topmost face over land, and the first wet face above the
# bathymetry — carry the identity row, which is how Υ(η) = Υ(-H) = 0 enters. `Υ` arrives holding
# Υᴳᴹ, is overwritten with the modified right-hand side on the way up, then with the solution on
# the way down; `t` holds the modified upper diagonal.
@inline function solve_transport_column!(Υ, t, i, j, grid, buoyancy, tracers,
                                         ℓx, ℓy, N², Δzᶜ, Δzᶠ, M, c★, N²₀)

    Nz = size(grid, 3)
    c² = squared_wave_speed(i, j, grid, buoyancy, tracers, ℓx, ℓy, N², Δzᶜ, Δzᶠ, M, c★, N²₀)

    @inbounds begin
        # Bottom face: always the identity row, so the sweep starts from a diagonal of one.
        t[i, j, 1] = zero(grid)
        Υ[i, j, 1] = zero(grid)

        for k in 2:Nz+1
            interior = !peripheral_node(i, j, k, grid, ℓx, ℓy, Face())

            N²ᵏ = max(finite_or_zero(N²(i, j, k, grid, buoyancy, tracers), grid), N²₀)
            aₖ  = c² / (Δzᶜ(i, j, k-1, grid) * Δzᶠ(i, j, k, grid))
            cₖ  = c² / (Δzᶜ(i, j, k,   grid) * Δzᶠ(i, j, k, grid))
            bₖ  = - (aₖ + cₖ) - N²ᵏ
            dₖ  = - N²ᵏ * Υ[i, j, k]

            # Identity row on the boundary of the fluid.
            aₖ = ifelse(interior, aₖ, zero(grid))
            cₖ = ifelse(interior, cₖ, zero(grid))
            bₖ = ifelse(interior, bₖ, one(grid))
            dₖ = ifelse(interior, dₖ, zero(grid))

            β = bₖ - aₖ * t[i, j, k-1]
            t[i, j, k] = cₖ / β
            Υ[i, j, k] = (dₖ - aₖ * Υ[i, j, k-1]) / β
        end

        for k in Nz:-1:1
            Υ[i, j, k] -= t[i, j, k] * Υ[i, j, k+1]
        end
    end

    return nothing
end

@kernel function _solve_transport_boundary_value_problem!(Υˣ, Υʸ, t, grid, buoyancy, tracers, M, c★, N²₀)
    i, j = @index(Global, NTuple)

    solve_transport_column!(Υˣ, t, i, j, grid, buoyancy, tracers,
                            Face(), Center(), N²ᶠᶜᶠ, Δzᶠᶜᶜ, Δzᶠᶜᶠ, M, c★, N²₀)

    solve_transport_column!(Υʸ, t, i, j, grid, buoyancy, tracers,
                            Center(), Face(), N²ᶜᶠᶠ, Δzᶜᶠᶜ, Δzᶜᶠᶠ, M, c★, N²₀)
end

#####
##### Eddy-induced velocity, (u★, w★) = (∂z Υ, -∇ ⋅ Υ)
#####

@kernel function _compute_bvp_eddy_velocities!(uₑ, vₑ, wₑ, grid, Υˣ, Υʸ)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        uₑ[i, j, k] = δzᵃᵃᶜ(i, j, k, grid, Υˣ) * Δz⁻¹ᶠᶜᶜ(i, j, k, grid)
        vₑ[i, j, k] = δzᵃᵃᶜ(i, j, k, grid, Υʸ) * Δz⁻¹ᶜᶠᶜ(i, j, k, grid)

        wˣ = δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶠ, Υˣ)
        wʸ = δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶠ, Υʸ)

        wₑ[i, j, k] = - (wˣ + wʸ) * Az⁻¹ᶜᶜᶠ(i, j, k, grid)
    end
end
