using Oceananigans.OutputReaders: FieldTimeSeries, InMemory
using Oceananigans.Fields: interior
using Oceananigans.Operators: Δxᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶜᶠᶜ, Δzᶠᶜᶜ

"""
    StraitSection

A rectangular strait section on a tripolar/ORCA grid.

* `i`, `j` are 1-based index ranges into the interior grid.
* `axis` is `:v` for a zonal (constant-`j`) section, where transport is
  ``∑ vₒ \\, Δx \\, Δz``, and `:u` for a meridional (constant-`i`) section,
  where transport is ``∑ uₒ \\, Δy \\, Δz``.
"""
struct StraitSection
    i :: UnitRange{Int}
    j :: UnitRange{Int}
    axis :: Symbol
end

# Per-configuration section indices. The half-degree indices are derived
# from a 720x360 TripolarGrid; the ORCA indices from ORCAGrid(ORCA1()).
# Bering, Drake and ITF are picked at the cells closest to standard
# observational sections (Bering Strait ~66°N/169°W, Drake ~67°W/57°S,
# ITF ~110°-130°E/8.5°S).
strait_sections(::Val{:halfdegree}) = (
    bering = StraitSection(212:218, 314:314, :v),
    drake  = StraitSection(447:447,  32:54,  :u),
    itf    = StraitSection( 83:122, 154:154, :v),
)

strait_sections(::Val{:orca}) = (
    bering = StraitSection(112:118, 251:251, :v),
    drake  = StraitSection(221:221,  53:71,  :u),
    itf    = StraitSection( 39:58,  130:130, :v),
)

strait_sections(config::Symbol) = strait_sections(Val(config))

"""
    strait_transports(config::Symbol, fields_file::AbstractString;
                      backend = InMemory(10),
                      start_time = 0, stop_time = Inf)

Compute time series of volume transport (Sv) through the Bering Strait,
Drake Passage and the Indonesian Throughflow from the offline 3-D output
`fields_file` (typically `<prefix>_fields.jld2`).

Dispatches on `config`: `:halfdegree` for the 720x360 TripolarGrid,
`:orca` for the ORCA1 mesh.

Returns `(; bering, drake, itf, time)` where each transport is a
`Vector{Float64}` in Sverdrups.
"""
function strait_transports(config::Symbol, fields_file::AbstractString;
                           backend = InMemory(10),
                           start_time = 0,
                           stop_time = Inf)

    sections = strait_sections(config)

    u_fts = FieldTimeSeries(fields_file, "uo"; backend = deepcopy(backend))
    v_fts = FieldTimeSeries(fields_file, "vo"; backend = deepcopy(backend))
    grid  = u_fts.grid

    times = collect(u_fts.times)
    Nt = length(times)
    bering = zeros(Nt)
    drake  = zeros(Nt)
    itf    = zeros(Nt)

    for n in 1:Nt
        u_int = interior(u_fts[n])
        v_int = interior(v_fts[n])
        bering[n] = section_volume_flux(grid, u_int, v_int, sections.bering) * 1e-6
        drake[n]  = section_volume_flux(grid, u_int, v_int, sections.drake)  * 1e-6
        itf[n]    = section_volume_flux(grid, u_int, v_int, sections.itf)    * 1e-6
    end

    in_window = (times .>= start_time) .& (times .<= stop_time)
    return (bering = bering[in_window],
            drake  = drake[in_window],
            itf    = itf[in_window],
            time   = times[in_window])
end

function section_volume_flux(grid, u_int, v_int, section::StraitSection)
    Nz = size(u_int, 3)
    total = 0.0

    if section.axis == :v
        for j in section.j, i in section.i, k in 1:Nz
            Δx = Δxᶜᶠᶜ(i, j, k, grid)
            Δz = Δzᶜᶠᶜ(i, j, k, grid)
            total += v_int[i, j, k] * Δx * Δz
        end
    elseif section.axis == :u
        for j in section.j, i in section.i, k in 1:Nz
            Δy = Δyᶠᶜᶜ(i, j, k, grid)
            Δz = Δzᶠᶜᶜ(i, j, k, grid)
            total += u_int[i, j, k] * Δy * Δz
        end
    else
        throw(ArgumentError("section.axis must be :u or :v, got $(section.axis)"))
    end

    return total
end
