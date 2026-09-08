using Adapt: Adapt
using Oceananigans.TurbulenceClosures: TurbulenceClosures

# The triad closure evaluates one tapering factor per cell, from a slope that is not the slope it
# multiplies: `tapering_factorᶜᶜᶜ` divides an interpolated `∂x b` by an x-averaged `∂z b`, while
# each triad divides a one-sided `∂x b` by a single raw `∂z b`. Where one corner of the stencil is
# weakly stratified and its neighbour is not, the averaged denominator stays healthy, the factor
# stays at 1, and the triad slope runs away. Measured on the ORCA initial condition, `ϵκR₃₃`
# reaches 7.3e4 m² s⁻¹ against the 2 κ Sₘ² = 0.16 the limiter is meant to guarantee, and the
# explicit R₃₁/R₃₂ terms that carry the same slopes blow the run up within a handful of steps.
# `IsopycnalSkewSymmetricDiffusivity` is immune because its tapering factor and its rotation
# tensor are built from the same derivatives.
#
# `TriadSlopeTapering` tapers each triad on its own slope, which restores `ϵ |S| ≤ Sₘ` and
# `ϵ S² ≤ Sₘ²` triad by triad. The same scalar multiplies a triad's `∂z c` term in the horizontal
# flux and its `∂x c` term in the vertical flux, so the skew tensor stays antisymmetric.
#
# TODO: this belongs in Oceananigans, in `isopycnal_skew_symmetric_diffusivity_with_triads.jl` —
# `ϵx⁺⁺` … `ϵy⁻⁻` should taper on their own slope rather than call `tapering_factorᶜᶜᶜ`.
# Upstreamed in CliMA/Oceananigans.jl#5938.

"""
    struct TriadSlopeTapering{L}

Slope limiter for `TriadIsopycnalSkewSymmetricDiffusivity` that applies the wrapped `limiter` to
each triad slope separately, rather than once per cell to an interpolated slope.
"""
struct TriadSlopeTapering{L}
    limiter :: L
end

Adapt.adapt_structure(to, tapering::TriadSlopeTapering) = TriadSlopeTapering(Adapt.adapt(to, tapering.limiter))

@inline TurbulenceClosures.tapering_factor(Sx, Sy, tapering::TriadSlopeTapering) =
    TurbulenceClosures.tapering_factor(Sx, Sy, tapering.limiter)

@inline triad_taper(S, sl::TriadSlopeTapering) = TurbulenceClosures.tapering_factor(S, zero(S), sl.limiter)

@inline TurbulenceClosures.ϵx⁺⁺(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_x(i+1, i, j, k, k+1, grid) * triad_taper(TurbulenceClosures.Sx⁺⁺(i, j, k, grid, b, C), sl)
@inline TurbulenceClosures.ϵx⁺⁻(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_x(i+1, i, j, k, k,   grid) * triad_taper(TurbulenceClosures.Sx⁺⁻(i, j, k, grid, b, C), sl)
@inline TurbulenceClosures.ϵx⁻⁺(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_x(i,   i, j, k, k+1, grid) * triad_taper(TurbulenceClosures.Sx⁻⁺(i, j, k, grid, b, C), sl)
@inline TurbulenceClosures.ϵx⁻⁻(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_x(i,   i, j, k, k,   grid) * triad_taper(TurbulenceClosures.Sx⁻⁻(i, j, k, grid, b, C), sl)

@inline TurbulenceClosures.ϵy⁺⁺(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_y(i, j+1, j, k, k+1, grid) * triad_taper(TurbulenceClosures.Sy⁺⁺(i, j, k, grid, b, C), sl)
@inline TurbulenceClosures.ϵy⁺⁻(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_y(i, j+1, j, k, k,   grid) * triad_taper(TurbulenceClosures.Sy⁺⁻(i, j, k, grid, b, C), sl)
@inline TurbulenceClosures.ϵy⁻⁺(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_y(i, j,   j, k, k+1, grid) * triad_taper(TurbulenceClosures.Sy⁻⁺(i, j, k, grid, b, C), sl)
@inline TurbulenceClosures.ϵy⁻⁻(i, j, k, grid, sl::TriadSlopeTapering, b, C) =
    TurbulenceClosures.triad_mask_y(i, j,   j, k, k,   grid) * triad_taper(TurbulenceClosures.Sy⁻⁻(i, j, k, grid, b, C), sl)
