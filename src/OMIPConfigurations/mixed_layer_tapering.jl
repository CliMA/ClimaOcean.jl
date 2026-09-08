using Adapt: Adapt
using KernelAbstractions: @index, @kernel
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.BuoyancyFormulations: ∂xᵣ_b, ∂yᵣ_b, ∂z_b, buoyancy_perturbationᶜᶜᶜ
using Oceananigans.Grids: Center, Face, inactive_node, znode
using Oceananigans.Operators: Δzᶜᶜᶜ, ∂x_zᶠᶜᶜ, ∂y_zᶠᶜᶜ, ∂x_zᶜᶠᶜ, ∂y_zᶜᶠᶜ, ∂x_zᶜᶜᶠ, ∂y_zᶜᶜᶠ,
                              ℑxyᶠᶜᵃ, ℑxzᶠᵃᶜ, ℑxyᶜᶠᵃ, ℑyzᵃᶠᶜ, ℑxzᶜᵃᶠ, ℑyzᵃᶜᶠ
using Oceananigans.TurbulenceClosures: TurbulenceClosures, IsopycnalSkewSymmetricDiffusivity, calc_tapering
using Oceananigans.Utils: KernelParameters, launch!

# Mixed-layer tapering of the isopycnal slopes (Danabasoglu, Ferrari & McWilliams 2008;
# NEMO's `ldfslp`): on top of the Gerdes et al. (1991) slope clip, the eddy fluxes are ramped
# linearly to zero from the local mixed-layer base to the surface. Without this, the GM
# overturning acts at full strength through weakly stratified winter columns — the clipped
# slope sits at the cap through the whole column — and restratifies deep-convection sites far
# more efficiently than in models that taper (which is why NEMO survives the same Treguier
# coefficient that kills the Labrador Sea here). The mixed-layer depth is a plain field
# refreshed each step by `RefreshMixedLayerTapering` with the dBM (2004) buoyancy criterion,
# mirroring the coefficient callbacks (and sharing their one-step pickup transient).

"""
    struct MixedLayerTapering{FT, F}

Slope limiter for `IsopycnalSkewSymmetricDiffusivity` combining the Gerdes-style `max_slope`
clip with a linear ramp of the tapering factor from the local mixed-layer base to the surface.
"""
struct MixedLayerTapering{FT, F}
    max_slope :: FT                   # duck-types FluxTapering for `calc_tapering`
    mixed_layer_depth :: F            # (Center, Center, Nothing), metres, positive down
end

function MixedLayerTapering(grid; max_slope = 1e-2)
    FT = eltype(grid)
    return MixedLayerTapering(convert(FT, max_slope), Field{Center, Center, Nothing}(grid))
end

Adapt.adapt_structure(to, limiter::MixedLayerTapering) =
    MixedLayerTapering(limiter.max_slope, Adapt.adapt(to, limiter.mixed_layer_depth))

# 0 at the surface, 1 at and below the mixed-layer base
@inline function mixed_layer_ramp(i, j, k, grid, limiter)
    z = znode(i, j, k, grid, Center(), Center(), Face())
    h = @inbounds limiter.mixed_layer_depth[i, j, 1]
    return clamp(-z / max(h, convert(eltype(grid), 1)), zero(grid), one(grid))
end

const MixedLayerTaperedISSD =
    IsopycnalSkewSymmetricDiffusivity{<:Any, <:Any, <:Any, <:Any, <:Any, <:MixedLayerTapering}

# The three position-aware tapering factors of the diffusive fluxes (and the implicit R₃₃
# precomputation, which reuses the ᶜᶜᶠ one): the Gerdes factor from `calc_tapering` times the
# mixed-layer ramp. Bodies mirror the untapered originals in
# `isopycnal_skew_symmetric_diffusivity.jl`.
@inline function TurbulenceClosures.tapering_factorᶠᶜᶜ(i, j, k, grid, closure::MixedLayerTaperedISSD, tracers, buoyancy)
    by   = ℑxyᶠᶜᵃ(i, j, k, grid, ∂yᵣ_b, buoyancy, tracers)
    bz   = ℑxzᶠᵃᶜ(i, j, k, grid, ∂z_b,  buoyancy, tracers)
    bx   =  ∂xᵣ_b(i, j, k, grid, buoyancy, tracers)
    ∂x_z = ∂x_zᶠᶜᶜ(i, j, k, grid)
    ∂y_z = ∂y_zᶠᶜᶜ(i, j, k, grid)
    ϵ = calc_tapering(bx, by, bz, ∂x_z, ∂y_z, grid, closure.isopycnal_tensor, closure.slope_limiter)
    return ϵ * mixed_layer_ramp(i, j, k, grid, closure.slope_limiter)
end

@inline function TurbulenceClosures.tapering_factorᶜᶠᶜ(i, j, k, grid, closure::MixedLayerTaperedISSD, tracers, buoyancy)
    bx   = ℑxyᶜᶠᵃ(i, j, k, grid, ∂xᵣ_b, buoyancy, tracers)
    bz   = ℑyzᵃᶠᶜ(i, j, k, grid, ∂z_b,  buoyancy, tracers)
    by   =  ∂yᵣ_b(i, j, k, grid, buoyancy, tracers)
    ∂x_z = ∂x_zᶜᶠᶜ(i, j, k, grid)
    ∂y_z = ∂y_zᶜᶠᶜ(i, j, k, grid)
    ϵ = calc_tapering(bx, by, bz, ∂x_z, ∂y_z, grid, closure.isopycnal_tensor, closure.slope_limiter)
    return ϵ * mixed_layer_ramp(i, j, k, grid, closure.slope_limiter)
end

@inline function TurbulenceClosures.tapering_factorᶜᶜᶠ(i, j, k, grid, closure::MixedLayerTaperedISSD, tracers, buoyancy)
    bx   = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᵣ_b, buoyancy, tracers)
    by   = ℑyzᵃᶜᶠ(i, j, k, grid, ∂yᵣ_b, buoyancy, tracers)
    bz   =  ∂z_b(i, j, k, grid, buoyancy, tracers)
    ∂x_z = ∂x_zᶜᶜᶠ(i, j, k, grid)
    ∂y_z = ∂y_zᶜᶜᶠ(i, j, k, grid)
    ϵ = calc_tapering(bx, by, bz, ∂x_z, ∂y_z, grid, closure.isopycnal_tensor, closure.slope_limiter)
    return ϵ * mixed_layer_ramp(i, j, k, grid, closure.slope_limiter)
end

# Advective-formulation slope functions: same ramp on top of the Gerdes magnitude clip.
@inline function TurbulenceClosures.ϵSxᶠᶜᶠ(i, j, k, grid, limiter::MixedLayerTapering, b, C)
    Sx = TurbulenceClosures.Sxᶠᶜᶠ(i, j, k, grid, b, C)
    ϵ  = TurbulenceClosures.tapering_factor(Sx, zero(grid), limiter)
    return ϵ * mixed_layer_ramp(i, j, k, grid, limiter) * Sx
end

@inline function TurbulenceClosures.ϵSyᶜᶠᶠ(i, j, k, grid, limiter::MixedLayerTapering, b, C)
    Sy = TurbulenceClosures.Syᶜᶠᶠ(i, j, k, grid, b, C)
    ϵ  = TurbulenceClosures.tapering_factor(zero(grid), Sy, limiter)
    return ϵ * mixed_layer_ramp(i, j, k, grid, limiter) * Sy
end

# Position-free magnitude clip, duck-typing FluxTapering
@inline TurbulenceClosures.tapering_factor(Sx, Sy, limiter::MixedLayerTapering) =
    min(one(Sx), limiter.max_slope^2 / (Sx^2 + Sy^2 + convert(typeof(Sx), 1e-40)))

#####
##### Mixed-layer depth refresh (dBM 2004 buoyancy criterion, Δb = 2.87e-4 m/s²)
#####

@kernel function _compute_tapering_mixed_layer_depth!(h, grid, buoyancy, tracers, Δb)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)
    bₛ = buoyancy_perturbationᶜᶜᶜ(i, j, Nz, grid, buoyancy, tracers)
    depth = zero(grid)
    found = zero(grid)
    for k in Nz:-1:1
        inactive = inactive_node(i, j, k, grid, Center(), Center(), Center())
        b = buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, buoyancy, tracers)
        z = znode(i, j, k, grid, Center(), Center(), Center())
        below = ifelse((bₛ - b > Δb) & (found == 0) & !inactive, one(grid), zero(grid))
        depth = ifelse(below == 1, -z, depth)
        found = max(found, below)
        depth = ifelse((found == 0) & !inactive, -z, depth)   # track column depth as fallback
    end
    @inbounds h[i, j, 1] = depth
end

"""
    compute_tapering_mixed_layer_depth!(limiter, ocean_model; Δb = 2.87e-4)

Refresh the limiter's mixed-layer-depth field from the current buoyancy field.
"""
function compute_tapering_mixed_layer_depth!(limiter::MixedLayerTapering, ocean_model; Δb = 2.87e-4)
    grid = ocean_model.grid
    launch!(architecture(grid), grid, :xy, _compute_tapering_mixed_layer_depth!,
            limiter.mixed_layer_depth, grid, ocean_model.buoyancy.formulation,
            ocean_model.tracers, convert(eltype(grid), Δb))

    # The ramp makes the taper factor depend on the depth at a single column, so the two stencil
    # representations of a tripolar fold face see different factors and the closure leaks salt
    # (−0.33 Gt/yr in production — 200× below the old coefficient leak, but nonzero). A zonally
    # uniform depth across the fold band restores single-valuedness, same argument as
    # `match_fold_rows!` for the coefficients.
    homogenize_fold_band!(limiter.mixed_layer_depth, grid)

    fill_halo_regions!(limiter.mixed_layer_depth)
    return nothing
end

homogenize_fold_band!(depth_field, grid; rows = 5) =
    launch!(architecture(grid), grid, KernelParameters((1,), (0,)),
            _homogenize_fold_band!, depth_field, grid, rows)

@kernel function _homogenize_fold_band!(h, grid, rows)
    @index(Global)
    Nx = size(grid, 1)
    Ny = size(grid, 2)
    μ = zero(grid)
    for i in 1:Nx
        @inbounds μ += h[i, Ny - rows, 1]
    end
    μ = μ / Nx
    for j in Ny-rows+1:Ny, i in 1:Nx
        @inbounds h[i, j, 1] = μ
    end
end

struct RefreshMixedLayerTapering{L}
    limiter :: L
end

(r::RefreshMixedLayerTapering)(sim) = compute_tapering_mixed_layer_depth!(r.limiter, sim.model.ocean.model)
