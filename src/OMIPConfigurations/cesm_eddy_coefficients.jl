using KernelAbstractions: @index, @kernel
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.BuoyancyFormulations: ∂z_b
using Oceananigans.Grids: Center, Face, inactive_node, znode
using Oceananigans.Utils: launch!

# CESM/POP's stratification-dependent vertical structure for the GM and Redi coefficients
# (Danabasoglu & Marshall 2007): κ(z) = κ_ref · clamp(N²(z)/N²_ref, r_min, 1), where N²_ref is
# the reference stratification just below the surface diabatic layer — approximated here as the
# maximum N² over the upper `reference_depth` of the column. The coefficient therefore stays at
# its reference value in the stratified thermocline and collapses to `r_min · κ_ref` both below
# it and throughout weakly stratified (convecting) columns, which is what lets wintertime deep
# convection survive the closure while the subtropical pycnocline keeps its full restraint.
#
# Complementary to the Treguier option in `nemo_eddy_coefficients.jl`: that one is horizontally
# selective but depth-uniform; this one carries the vertical structure a constant coefficient
# cannot. Both plug the same (Center, Center, Face) coefficient fields into
# `IsopycnalSkewSymmetricDiffusivity` and are refreshed each time step by a callback.
#
# Pickup caveat (shared with the NEMO coefficients): the fields are primed at build time from
# the pre-pickup initial state, and `run!` only invokes time-step callbacks up front when
# `clock.iteration == 0`, so the first step after a checkpoint pickup uses coefficients one
# state behind. A one-step transient.

"""
    struct CESMEddyCoefficients{F, FT}

Stratification-dependent GM and Redi coefficients following CESM's Danabasoglu & Marshall
(2007) vertical structure. Both fields sit at the `IsopycnalSkewSymmetricDiffusivity`
coefficient location so the closure reads them without interpolation, and both are recomputed
in place from the model state, so stepping allocates nothing.
"""
struct CESMEddyCoefficients{F, FT}
    skew_coefficient      :: F
    symmetric_coefficient :: F
    parameters            :: NTuple{5, FT}   # κˢ₀, κʳ₀, minimum fraction, reference depth, minimum N²_ref
end

function CESMEddyCoefficients(grid;
                              reference_skew_coefficient = 1000,
                              reference_symmetric_coefficient = 1000,
                              minimum_fraction = 0.1,               # D&M (2007) lower clip
                              reference_depth = 500,                # m, window for N²_ref
                              minimum_reference_stratification = 1e-9)  # s⁻²; below ⇒ unstratified column

    FT = eltype(grid)
    parameters = (convert(FT, reference_skew_coefficient),
                  convert(FT, reference_symmetric_coefficient),
                  convert(FT, minimum_fraction),
                  convert(FT, reference_depth),
                  convert(FT, minimum_reference_stratification))

    return CESMEddyCoefficients(Field{Center, Center, Face}(grid),
                                Field{Center, Center, Face}(grid),
                                parameters)
end

@kernel function _compute_cesm_eddy_coefficients!(κˢ, κʳ, grid, buoyancy, tracers, parameters)
    i, j = @index(Global, NTuple)
    κˢ₀, κʳ₀, rmin, dref, N²min = parameters
    Nz = size(grid, 3)

    # Pass 1: N²_ref = maximum stratification over the interior faces of the upper `dref` meters.
    N²ref = zero(grid)
    wet = zero(grid)
    for k in 2:Nz
        inactive = inactive_node(i, j, k, grid, Center(), Center(), Face())
        z = znode(i, j, k, grid, Center(), Center(), Face())
        N² = max(0, finite_or_zero(∂z_b(i, j, k, grid, buoyancy, tracers), grid))
        in_window = ifelse(inactive, zero(grid), ifelse(z ≥ -dref, one(grid), zero(grid)))
        N²ref = max(N²ref, in_window * N²)
        wet = max(wet, in_window)
    end

    unstratified = N²ref < N²min

    # Pass 2: fill the interior faces; boundary faces copy their interior neighbor below/above.
    for k in 2:Nz
        N² = max(0, finite_or_zero(∂z_b(i, j, k, grid, buoyancy, tracers), grid))
        ratio = ifelse(unstratified, rmin, clamp(N² / max(N²ref, N²min), rmin, one(grid)))
        ratio = wet * ratio
        @inbounds κˢ[i, j, k] = κˢ₀ * ratio
        @inbounds κʳ[i, j, k] = κʳ₀ * ratio
    end

    @inbounds κˢ[i, j, 1] = κˢ[i, j, 2]
    @inbounds κʳ[i, j, 1] = κʳ[i, j, 2]
    @inbounds κˢ[i, j, Nz+1] = κˢ[i, j, Nz]
    @inbounds κʳ[i, j, Nz+1] = κʳ[i, j, Nz]
end

"""
    compute_cesm_eddy_coefficients!(coefficients, ocean_model)

Refresh both coefficient fields from the current buoyancy field.
"""
function compute_cesm_eddy_coefficients!(coefficients::CESMEddyCoefficients, ocean_model)
    grid = ocean_model.grid

    launch!(architecture(grid), grid, :xy, _compute_cesm_eddy_coefficients!,
            coefficients.skew_coefficient, coefficients.symmetric_coefficient, grid,
            ocean_model.buoyancy, fields(ocean_model), coefficients.parameters)

    match_fold_rows!(coefficients.symmetric_coefficient, coefficients.skew_coefficient, grid)

    # Without this the halos stay zero, so the two cells sharing a face across the periodic seam or the
    # tripolar fold evaluate the isopycnal flux with different κ. The face flux is then not antisymmetric
    # between them, the divergence stops telescoping, and the closure leaks salt at a steady rate.
    fill_halo_regions!(coefficients.skew_coefficient)
    fill_halo_regions!(coefficients.symmetric_coefficient)

    return nothing
end

struct RefreshCESMEddyCoefficients{C}
    coefficients :: C
end

(r::RefreshCESMEddyCoefficients)(sim) = compute_cesm_eddy_coefficients!(r.coefficients, sim.model.ocean.model)

# `:cesm` selects the Danabasoglu & Marshall coefficient for either diffusivity; anything else
# passes through to `IsopycnalSkewSymmetricDiffusivity` unchanged.
uses_cesm_eddy_coefficients(κ_skew, κ_symmetric) = (κ_skew === :cesm) | (κ_symmetric === :cesm)

resolve_cesm_coefficient(κ, ::Nothing, field_name) = κ
resolve_cesm_coefficient(κ, coefficients, field_name) =
    κ === :cesm ? getproperty(coefficients, field_name) : κ
