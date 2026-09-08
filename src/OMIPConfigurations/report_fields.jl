using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Operators: ζ₃ᶠᶠᶜ, ℑxᶜᵃᵃ, ℑyᵃᶜᵃ
using Oceananigans.Architectures: child_architecture
using Oceananigans.Fields: interpolate!
using NumericalEarth.DataWrangling: WOAAnnual
using NumericalEarth.Diagnostics: MixedLayerDepthField
using WorldOceanAtlasTools

@inline function speedᶜᶜᶜ(i, j, k, grid, u, v)
    û = ℑxᶜᵃᵃ(i, j, k, grid, u)
    v̂ = ℑyᵃᶜᵃ(i, j, k, grid, v)
    return sqrt(û^2 + v̂^2)
end

"""
    compute_report_fields(ocean; dataset = WOAAnnual())

Compute a `NamedTuple` of diagnostic fields from the current state of
`ocean`. Returns surface-level slices and zonal averages, plus
differences against the WOA climatology specified by `dataset`.

Returned fields:
- `SST`, `SSS`: 2-D surface temperature and salinity
- `spd`: surface speed sqrt(u^2 + v^2)
- `ζ`: surface vertical vorticity
- `MLD`: mixed-layer depth
- `T̄`, `S̄`: zonally averaged temperature and salinity (latitude × depth)
- `δT`, `δS`: SST/SSS minus WOA climatology
- `φ`: latitude coordinates of the zonal averages
- `z`: depth coordinates of the zonal averages
"""
function compute_report_fields(ocean; dataset = WOAAnnual())
    grid = ocean.model.grid
    arch = child_architecture(grid)
    Nz = size(grid, 3)

    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S

    SST = Array(interior(T, :, :, Nz))
    SSS = Array(interior(S, :, :, Nz))

    spd_op = KernelFunctionOperation{Center, Center, Center}(speedᶜᶜᶜ, grid, u, v)
    spd_field = Field(spd_op; indices = (:, :, Nz))
    compute!(spd_field)
    spd = Array(interior(spd_field, :, :, 1))

    ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
    ζ_field = Field(ζ_op; indices = (:, :, Nz))
    compute!(ζ_field)
    ζ = Array(interior(ζ_field, :, :, 1))

    h = MixedLayerDepthField(ocean.model.buoyancy, grid, ocean.model.tracers)
    compute!(h)
    MLD = Array(interior(h, :, :, 1))

    δT, δS = compute_woa_bias(grid, arch, T, S, Nz, dataset)

    return (; SST, SSS, spd, ζ, MLD, δT, δS)
end

"""
    compute_woa_bias(grid, arch, T, S, Nz, dataset)

Return `(δT, δS)`, the surface temperature and salinity differences
between the current state and the WOA climatology specified by
`dataset` (default `WOAAnnual()`).
"""
function compute_woa_bias(grid, arch, T, S, Nz, dataset)
    Tʷ = Field(Metadatum(:temperature; dataset), arch)
    Sʷ = Field(Metadatum(:salinity;    dataset), arch)

    Tᵢ = CenterField(grid)
    Sᵢ = CenterField(grid)
    interpolate!(Tᵢ, Tʷ)
    interpolate!(Sᵢ, Sʷ)

    δT = Array(interior(T, :, :, Nz)) .- Array(interior(Tᵢ, :, :, Nz))
    δS = Array(interior(S, :, :, Nz)) .- Array(interior(Sᵢ, :, :, Nz))

    return δT, δS
end
