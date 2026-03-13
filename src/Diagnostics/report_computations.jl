"""
    compute_report_fields(ocean; dataset = WOAAnnual())

Compute diagnostic fields from the current state of `ocean`.
Returns a `NamedTuple` with:
- `SST`, `SSS`: surface temperature and salinity (2D arrays)
- `|u|`: surface speed √(u² + v²) (2D array)
- `ζ`: surface vertical vorticity ζ₃ (2D array)
- `h`: mixed layer depth (2D array)
- `T̄`, `S̄`: zonally averaged temperature and salinity (2D arrays, latitude × depth)
- `δT`, `δS`: SST/SSS minus WOA climatology (2D arrays)
- `φ`: latitude coordinates for zonal averages
- `z`: depth coordinates for zonal averages
"""
function compute_report_fields(ocean; dataset = WOAAnnual())
    grid = ocean.model.grid
    arch = child_architecture(grid)
    Nz = size(grid, 3)

    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S

    # Surface fields
    SST = Array(interior(T, :, :, Nz))
    SSS = Array(interior(S, :, :, Nz))

    # Surface speed
    uₛ = Field(u; indices = (:, :, Nz))
    vₛ = Field(v; indices = (:, :, Nz))
    c = @at (Center, Center, Nothing) sqrt(uₛ^2 + vₛ^2)
    compute!(Field(c))
    spd = Array(interior(Field(c), :, :, 1))

    # Surface vorticity
    ζ_op = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
    compute!(Field(ζ_op; indices = (:, :, Nz)))
    ζ = Array(interior(Field(ζ_op; indices = (:, :, Nz)), :, :, 1))

    # Mixed layer depth
    h = MixedLayerDepthField(ocean.model.buoyancy, grid, ocean.model.tracers)
    compute!(h)
    MLD = Array(interior(h, :, :, 1))

    # Zonal averages
    T̄, S̄, φ, z = compute_zonal_averages(grid, T, S)

    # WOA bias
    δT, δS = compute_woa_bias(grid, arch, T, S, Nz, dataset)

    return (; SST, SSS, spd, ζ, MLD, T̄, S̄, δT, δS, φ, z)
end

# --- Zonal average helpers ---

const LLGrid = Oceananigans.Grids.LatitudeLongitudeGrid
const IBG    = Oceananigans.ImmersedBoundaries.ImmersedBoundaryGrid

underlying(grid::IBG) = grid.underlying_grid
underlying(grid)      = grid

function compute_zonal_averages(grid, T, S)
    ug = underlying(grid)
    return compute_zonal_averages_on(ug, grid, T, S)
end

# Direct zonal average on LatitudeLongitudeGrid
function compute_zonal_averages_on(::LLGrid, grid, T, S)
    T̄ = compute!(Field(Average(T, dims = 1)))
    S̄ = compute!(Field(Average(S, dims = 1)))
    φ = Array(ynodes(T̄))
    z = Array(znodes(T̄))
    return dropdims(Array(interior(T̄)), dims = 1),
           dropdims(Array(interior(S̄)), dims = 1), φ, z
end

# For TripolarGrid (or any non-LatLon grid): regrid to lat-lon first,
# then compute the zonal average on the regular grid.
function compute_zonal_averages_on(ug, grid, T, S)
    Nx, Ny, Nz = size(grid)

    φλ_grid = LatitudeLongitudeGrid(CPU();
                                    size = (Nx, Ny, 1),
                                    longitude = (0, 360),
                                    latitude = (-90, 90),
                                    z = (0, 1))

    src = CenterField(grid; indices = (:, :, 1))
    dst = CenterField(φλ_grid)
    rgd = Regridder(dst, src)

    T̄ = zeros(Ny, Nz)
    S̄ = zeros(Ny, Nz)

    for k in 1:Nz
        Tₖ = Array(interior(T, :, :, k))
        regrid!(vec(interior(dst)), rgd, vec(Tₖ))
        T̄[:, k] = dropdims(mean(interior(dst), dims = 1), dims = 1)

        Sₖ = Array(interior(S, :, :, k))
        regrid!(vec(interior(dst)), rgd, vec(Sₖ))
        S̄[:, k] = dropdims(mean(interior(dst), dims = 1), dims = 1)
    end

    φ = Array(ynodes(dst))
    z = Array(znodes(T))

    return T̄, S̄, φ, z
end

# --- WOA bias helper ---

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
