using KernelAbstractions: @index, @kernel
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.BuoyancyFormulations: ∂z_b
using Oceananigans.Grids: Center, Face, inactive_node, φnode, znode
using Oceananigans.TurbulenceClosures: FluxTapering
using Oceananigans.Utils: KernelParameters, launch!

# Treguier × Danabasoglu-Marshall hybrid: the horizontal structure of the NEMO coefficient
# (Ro² × baroclinic growth rate, capped) multiplied by the CESM vertical shape
# clamp(N²/N²_ref, r_min, 1). The horizontal factor concentrates the mixing in the subtropics
# and ACC and starves the subpolar interiors; the vertical shape confines it to the stratified
# part of each column, so the weakly stratified interior of the subpolar dome — the multi-year
# preconditioning for deep convection — is protected in every season. This addresses the
# knemo failure mode: Ro = 0.4 ∫N dz/|f| breathes with the seasonal thermocline, inflating the
# depth-uniform Treguier coefficient ~4× in the summer subpolar gyre (Labrador ~160 in March,
# ~700 in August at year 25), which erodes the dome each restratification season and bleeds
# the AMOC over a decade. The default `maximum_rossby_radius` is also tightened to 20 km
# (NEMO's 40 km is what lets the summer subpolar inflate to the cap).

"""
    struct HybridEddyCoefficients{F, L, FT}

GM and Redi coefficients with Treguier et al. (1997) horizontal structure and a
Danabasoglu & Marshall (2007) stratification-dependent vertical shape. Both fields sit at the
`IsopycnalSkewSymmetricDiffusivity` coefficient location and are recomputed in place from the
model state, so stepping allocates nothing.
"""
struct HybridEddyCoefficients{F, L, FT}
    skew_coefficient      :: F
    symmetric_coefficient :: F
    slope_limiter         :: L
    parameters            :: NTuple{6, FT}
end

function HybridEddyCoefficients(grid;
                                maximum_skew_coefficient = 1000,
                                reference_symmetric_coefficient = 1000,
                                minimum_rossby_radius = 2e3,
                                maximum_rossby_radius = 20e3,
                                minimum_fraction = 0.1,
                                reference_stratification = 1e-5,
                                slope_limiter = FluxTapering(1e-2))

    FT = eltype(grid)
    parameters = (convert(FT, maximum_skew_coefficient),
                  convert(FT, reference_symmetric_coefficient),
                  convert(FT, minimum_rossby_radius),
                  convert(FT, maximum_rossby_radius),
                  convert(FT, minimum_fraction),
                  convert(FT, reference_stratification))

    return HybridEddyCoefficients(Field{Center, Center, Face}(grid),
                                  Field{Center, Center, Face}(grid),
                                  slope_limiter, parameters)
end

@kernel function _compute_hybrid_eddy_coefficients!(aei, aht, grid, limiter, buoyancy, tracers, Ω, parameters)
    i, j = @index(Global, NTuple)
    aei0, aht0, Romin, Romax, rmin, N²₀ = parameters
    Nz = size(grid, 3)

    # Column integrals: ∫N dz and ∫N²S² dz for the Treguier factor. The vertical shape uses a
    # constant reference N²₀ (see header note), applied pointwise in the second loop.
    zn = zero(grid)
    zah = zero(grid)
    zhw = zero(grid)
    for k in 1:Nz
        inactive = inactive_node(i, j, k, grid, Center(), Center(), Face())
        Δz = ifelse(inactive, zero(grid), Δzᶜᶜᶠ(i, j, k, grid))
        N² = max(finite_or_zero(∂z_b(i, j, k, grid, buoyancy, tracers), grid), zero(grid))
        N² = ifelse(inactive, zero(grid), N²)
        zn  += sqrt(N²) * Δz
        zah += N² * squared_isopycnal_slope(i, j, k, grid, limiter, buoyancy, tracers) * Δz
        zhw += Δz
    end

    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    f = 2 * Ω * sind(φ)
    f20 = 2 * Ω * sind(convert(eltype(grid), 20))

    Ro = clamp(convert(eltype(grid), 2//5) * zn / max(abs(f), convert(eltype(grid), 1e-10)), Romin, Romax)
    growth_rate = sqrt(zah / max(zhw, one(grid)))

    # Tropical factor: GM scaled by it (→ 0 at the equator), Redi by its complement (→ aht0).
    tropical = min(one(grid), abs(f) / f20)

    aeiw = min(tropical * Ro^2 * growth_rate, aei0)
    ahtmin = convert(eltype(grid), 1//5) * aht0
    ahtw = max(ahtmin, aeiw) + (one(grid) - tropical) * (aht0 - ahtmin)
    ahtw = min(ahtw, max(aei0, aht0))

    dry = zhw == zero(grid)
    aeiw = ifelse(dry, zero(grid), aeiw)
    ahtw = ifelse(dry, zero(grid), ahtw)
    for k in 2:Nz
        N² = max(0, finite_or_zero(∂z_b(i, j, k, grid, buoyancy, tracers), grid))
        shape = clamp(N² / N²₀, rmin, one(grid))
        @inbounds aei[i, j, k] = aeiw * shape
        @inbounds aht[i, j, k] = ahtw * shape
    end

    @inbounds aei[i, j, 1] = aei[i, j, 2]
    @inbounds aht[i, j, 1] = aht[i, j, 2]
    @inbounds aei[i, j, Nz+1] = aei[i, j, Nz]
    @inbounds aht[i, j, Nz+1] = aht[i, j, Nz]
end

"""
    compute_hybrid_eddy_coefficients!(coefficients, ocean_model)

Refresh both coefficient fields from the current buoyancy field.
"""
function compute_hybrid_eddy_coefficients!(coefficients::HybridEddyCoefficients, ocean_model)
    grid = ocean_model.grid
    Ω = convert(eltype(grid), Oceananigans.defaults.planet_rotation_rate)

    launch!(architecture(grid), grid, :xy, _compute_hybrid_eddy_coefficients!,
            coefficients.skew_coefficient, coefficients.symmetric_coefficient, grid,
            coefficients.slope_limiter, ocean_model.buoyancy, fields(ocean_model),
            Ω, coefficients.parameters)

    match_fold_rows!(coefficients.symmetric_coefficient, coefficients.skew_coefficient, grid)

    fill_halo_regions!(coefficients.skew_coefficient)
    fill_halo_regions!(coefficients.symmetric_coefficient)

    return nothing
end

struct RefreshHybridEddyCoefficients{C}
    coefficients :: C
end

(r::RefreshHybridEddyCoefficients)(sim) = compute_hybrid_eddy_coefficients!(r.coefficients, sim.model.ocean.model)

uses_hybrid_eddy_coefficients(κ_skew, κ_symmetric) = (κ_skew === :hybrid) | (κ_symmetric === :hybrid)

resolve_hybrid_coefficient(κ, ::Nothing, field_name) = κ
resolve_hybrid_coefficient(κ, coefficients, field_name) =
    κ === :hybrid ? getproperty(coefficients, field_name) : κ
