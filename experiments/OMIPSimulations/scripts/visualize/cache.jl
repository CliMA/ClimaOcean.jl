# visualize/cache.jl
#
# Per-case lazy field cache with optional disk persistence.
#
# Each derived diagnostic is keyed by a `Symbol` and computed on first
# access via the `LOADERS` registry below; subsequent accesses hit the
# in-memory cache. Loaders may pull other fields recursively
# (e.g. `:sst_bias` ← `:sst` + `:woa_temperature`), so each figure file
# just calls `get_field(case_cache, :sst_bias)` and the loader DAG
# handles the rest. Across figures, every field is loaded at most once
# per case.
#
# Loaders flagged with `disk_cached(...)` additionally persist their
# result under `<output_dir>/diag_cache/<prefix>/<sym>.jld2`. The cache
# is keyed by the snapshot count of every source `FieldTimeSeries` it
# depends on, plus the case averaging window — so appending new
# snapshots invalidates the entry automatically while no-op file
# touches do not.
#
# Requires `common.jl` to have been included first.

#####
##### CaseCache
#####

mutable struct CaseCache
    label      :: String
    prefix     :: String
    run_dir    :: String
    start_time :: Float64
    stop_time  :: Float64
    case       :: NamedTuple
    fields     :: Dict{Symbol, Any}
end

function CaseCache(case::NamedTuple)
    prefix  = case.prefix
    run_dir = run_dir_for(prefix)
    start_time, stop_time = resolve_case_window(case, run_dir, prefix)
    return CaseCache(case.label, prefix, run_dir,
                     start_time, stop_time,
                     case, Dict{Symbol, Any}())
end

"""
    resolve_case_window(case, run_dir, prefix) -> (start_time, stop_time)

Translate the case's averaging-window declaration into absolute
`(start_time, stop_time)` in seconds.

Two modes, mutually exclusive:

- **absolute**: any of `start_time` / `stop_time` are present. Missing
  fields default to `0.0` / `Inf`.
- **relative**: only `years_from_end` is present. The function reads
  the last snapshot time from the case's surface JLD2 stem and returns
  `(end_time - years_from_end * year, end_time + 1)`. The disk cache
  key uses these absolute values, so extending the run invalidates the
  cache automatically.

Raises if both modes are mixed or neither is given.
"""
function resolve_case_window(case::NamedTuple, run_dir::AbstractString, prefix::AbstractString)
    has_abs = haskey(case, :start_time) || haskey(case, :stop_time)
    has_rel = haskey(case, :years_from_end)

    if has_abs && has_rel
        error("Case '$(get(case, :label, prefix))': specify either `start_time`/`stop_time` or `years_from_end`, not both.")
    end

    if has_abs
        st = haskey(case, :start_time) ? Float64(case.start_time) : 0.0
        sp = haskey(case, :stop_time)  ? Float64(case.stop_time)  : Inf
        return st, sp
    end

    if has_rel
        surface_file = find_first_file(run_dir, prefix, "surface")
        end_time = try
            last_jld2_timeseries_time(surface_file)
        catch err
            error("Case '$(get(case, :label, prefix))': could not read last snapshot time from surface file — cannot derive window from `years_from_end`. ($(sprint(showerror, err)))")
        end
        window = Float64(case.years_from_end) * years
        return end_time - window, end_time + 1.0
    end

    error("Case '$(get(case, :label, prefix))': must specify either `start_time`/`stop_time` or `years_from_end`.")
end

# Sentinel used to memoize loaders that legitimately return `nothing`
# (e.g. dBM climatology when the file is missing). Without this a
# successful "no result" would re-trigger the loader on every access.
const CACHE_MISS = :__cache_miss__

"""
    get_field(case_cache, sym)

Return the value cached under `sym` for `case_cache`, computing it via
`LOADERS[sym]` on first access.
"""
function get_field(case_cache::CaseCache, sym::Symbol)
    val = get(case_cache.fields, sym, CACHE_MISS)
    val === CACHE_MISS || return val
    haskey(LOADERS, sym) || error("No loader registered for :$sym")
    val = LOADERS[sym](case_cache)
    case_cache.fields[sym] = val
    return val
end

#####
##### Disk-side caching
#####

"""
    DiskCacheKey(snapshot_counts, start_time, stop_time)

Validation key for a disk-cached derived field. A cached value is
reusable only if all three components match the current state of the
source `FieldTimeSeries` and the case averaging window.

`snapshot_counts` is an `NTuple{N, Int}`, one entry per source FTS the
loader depends on. Use `N = 0` for fields that depend only on the grid
or external climatologies (always-valid cache once written).
"""
struct DiskCacheKey{N}
    snapshot_counts :: NTuple{N, Int}
    start_time      :: Float64
    stop_time       :: Float64
end

Base.:(==)(a::DiskCacheKey{N}, b::DiskCacheKey{N}) where N =
    a.snapshot_counts == b.snapshot_counts &&
    a.start_time      == b.start_time      &&
    a.stop_time       == b.stop_time

# Cross-arity fallback. Required: keys whose `N` differs are not
# comparable via the parametric same-N method above (Julia would fall
# back to `===`, returning false only by accident). Without this method
# a `disk_cached` loader whose `source_fts_syms` arity changed across
# sessions would compare unequal *for the right reason* but with no
# helpful diagnostic — `explain_key_mismatch` also dispatches here.
Base.:(==)(::DiskCacheKey, ::DiskCacheKey) = false

"""
    diag_cache_dir(case_cache)

Per-case directory holding the JLD2 disk cache for postprocessed
diagnostics. Lives under `output_dir` so the run directory itself
stays read-only.
"""
diag_cache_dir(case_cache::CaseCache) = joinpath(output_dir, "diag_cache", case_cache.prefix)

diag_cache_path(case_cache::CaseCache, sym::Symbol) = joinpath(diag_cache_dir(case_cache), string(sym) * ".jld2")

# Maps every raw-FTS symbol (e.g. `:to_h_fts`) to the file symbol it
# lives in (`:averages_file`). Filled alongside each FTS loader
# registration below; consulted to validate disk caches via JLD2
# metadata only, never via a full `FieldTimeSeries` build.
const FTS_DISK_PATH_SYM = Dict{Symbol, Symbol}()

"""
    current_disk_cache_key(case_cache, fts_syms)

Validation key for a disk-cached derived field whose value depends on
the FTS named by each symbol in `fts_syms`. Reads `Nt` per source via
`total_jld2_timeseries_snapshot_count`, so no `FieldTimeSeries` is
constructed. Pass `()` for fields independent of model output.
"""
function current_disk_cache_key(c::CaseCache, fts_syms::Tuple{Vararg{Symbol}})
    Nts = ntuple(length(fts_syms)) do i
        fts_sym  = fts_syms[i]
        file_sym = get(FTS_DISK_PATH_SYM, fts_sym, nothing)
        isnothing(file_sym) && error("No JLD2 stem registered for FTS :$fts_sym")
        return total_jld2_timeseries_snapshot_count(get_field(c, file_sym))
    end
    return DiskCacheKey(Nts, c.start_time, c.stop_time)
end

current_disk_cache_key(c::CaseCache, fts_sym::Symbol) = current_disk_cache_key(c, (fts_sym,))

"""
    read_disk_cache(path)

Return `(value, key)` if `path` exists and parses cleanly; return
`nothing` otherwise. A failed read logs a warning and falls through
so the loader recomputes.
"""
function read_disk_cache(path)
    isfile(path) || return nothing
    try
        return JLD2.load(path, "value", "key")
    catch err
        @warn "  Disk cache read failed; will recompute" path=path error=sprint(showerror, err)
        return nothing
    end
end

"""
    write_disk_cache(path, value, key)

Persist `value` and its validation `key` to `path` (creating
directories as needed). Returns `value` so the caller can chain.
"""
function write_disk_cache(path, value, key)
    mkpath(dirname(path))
    try
        JLD2.jldsave(path; value, key)
    catch err
        @warn "  Disk cache write failed; continuing without persistence" path=path error=sprint(showerror, err)
    end
    return value
end

"""
    explain_key_mismatch(stored, current)

Return a short string describing what changed between `stored` and
`current` validation keys. Used to log *why* a disk cache entry is
being invalidated.
"""
function explain_key_mismatch(stored::DiskCacheKey{N}, current::DiskCacheKey{N}) where N
    diffs = String[]
    stored.snapshot_counts == current.snapshot_counts ||
        push!(diffs, "snapshots $(stored.snapshot_counts) → $(current.snapshot_counts)")
    stored.start_time == current.start_time ||
        push!(diffs, "start_time $(stored.start_time) → $(current.start_time)")
    stored.stop_time == current.stop_time ||
        push!(diffs, "stop_time $(stored.stop_time) → $(current.stop_time)")
    return join(diffs, ", ")
end

# Different N → keys aren't even comparable; report it explicitly.
explain_key_mismatch(::DiskCacheKey, ::DiskCacheKey) = "source-FTS arity changed"

"""
    disk_cached(loader, sym; source_fts_syms = ())

Wrap `loader(case_cache)` so its result is persisted under the
per-case JLD2 cache file for `sym`. Subsequent loads in any session
reuse the cached value when the validation key matches; otherwise
the loader is rerun and the cache is overwritten.

`source_fts_syms` is either a single FTS symbol or a tuple of them.
Use `()` for fields that depend only on the grid or climatologies.
"""
function disk_cached(loader::Function, sym::Symbol;
                     source_fts_syms::Union{Symbol, Tuple{Vararg{Symbol}}} = ())
    sources = source_fts_syms isa Symbol ? (source_fts_syms,) : source_fts_syms
    return function (c::CaseCache)
        path    = diag_cache_path(c, sym)
        new_key = current_disk_cache_key(c, sources)
        cached  = read_disk_cache(path)
        if !isnothing(cached)
            value, stored_key = cached
            if stored_key == new_key
                @info "  $(c.label): :$sym ← disk cache"
                return value
            end
            reason = stored_key isa DiskCacheKey ?
                     explain_key_mismatch(stored_key, new_key) :
                     "corrupt or unrecognized key"
            @info "  $(c.label): :$sym ← invalidated ($reason)"
        end
        @info "  $(c.label): :$sym ← computing"
        return write_disk_cache(path, loader(c), new_key)
    end
end

#####
##### Loader registry
#####
#
# Each entry is a function `(case_cache::CaseCache) -> value`. The
# returned value is cached verbatim, so loaders that return composite
# results (e.g. monthly bins tuple) can be unpacked by downstream
# loaders.

const LOADERS = Dict{Symbol, Function}()

#####
##### File paths
#####

LOADERS[:surface_file]  = c -> find_first_file(c.run_dir, c.prefix, "surface")
LOADERS[:fields_file]   = c -> find_first_file(c.run_dir, c.prefix, "fields")
LOADERS[:averages_file] = c -> find_first_file(c.run_dir, c.prefix, "averages")

#####
##### Raw FieldTimeSeries
#####

# Every raw `:..._fts` symbol → (file symbol, JLD2 variable name).
# One source of truth; consulted to register both the loader and the
# disk-cache file lookup.
const FTS_VARS = (
    surface_file  = ((:tos_fts,    "tos"),    (:sos_fts,   "sos"),
                     (:zos_fts,    "zos"),    (:mld_fts,   "mlotst"),
                     (:hfds_fts,   "hfds"),   (:wfo_fts,   "wfo"),
                     (:sic_fts,    "siconc"), (:zossq_fts, "zossq"),
                     (:tauuo_fts,  "tauuo"),  (:tauvo_fts, "tauvo"),
                     (:sithick_fts, "sithick")),
    fields_file   = ((:to_fts, "to"), (:so_fts, "so"), (:bo_fts, "bo"),
                     (:uo_fts, "uo"), (:vo_fts, "vo"),
                     (:uvol_fts, "uvol"), (:vvol_fts, "vvol")),
    averages_file = ((:tosga_fts, "tosga"), (:soga_fts, "soga"),
                     (:to_h_fts,  "to_h"),  (:so_h_fts,  "so_h")),
)

for (file_sym, mappings) in pairs(FTS_VARS), (sym, var) in mappings
    FTS_DISK_PATH_SYM[sym] = file_sym
    LOADERS[sym] = let v = var, f = file_sym
        c -> FieldTimeSeries(get_field(c, f), v; backend = deepcopy(FTS_BACKEND))
    end
end

#####
##### Grid + masks
#####

LOADERS[:grid]          = c -> total_jld2_serialized_grid(get_field(c, :surface_file), "tos")
LOADERS[:Nx]            = c -> size(get_field(c, :grid), 1)
LOADERS[:Ny]            = c -> size(get_field(c, :grid), 2)
LOADERS[:Nz]            = c -> size(get_field(c, :grid), 3)
LOADERS[:land]          = c -> build_land_mask(get_field(c, :grid))
LOADERS[:ocean_mask_3d] = c -> build_ocean_mask_3d(get_field(c, :grid))
LOADERS[:depth]         = c -> collect(znodes(get_field(c, :grid), Center()))

#####
##### Surface time means (with land masking)
#####

function time_mean_drop(c, fts_sym)
    fts = get_field(c, fts_sym)
    return dropdims(compute_time_mean(fts; start_time = c.start_time, stop_time = c.stop_time); dims = 3)
end

# Apply land mask in place; `nothing` is passed through (for optional
# climatology fields). The `maybe_*` helpers compose this with a few
# common wrap-or-skip patterns to remove `isnothing(...) ? nothing : …`
# boilerplate from the loader registry below.
masked!(c, ::Nothing) = nothing
masked!(c, value)     = (mask_land!(value, get_field(c, :land)); value)

maybe_diff!(c, a_sym, b_sym) = let b = get_field(c, b_sym)
    isnothing(b) ? nothing : masked!(c, get_field(c, a_sym) .- b)
end

LOADERS[:sst] = disk_cached(:sst; source_fts_syms = :tos_fts) do c
    masked!(c, time_mean_drop(c, :tos_fts))
end
LOADERS[:sss] = disk_cached(:sss; source_fts_syms = :sos_fts) do c
    masked!(c, time_mean_drop(c, :sos_fts))
end
LOADERS[:ssh] = disk_cached(:ssh; source_fts_syms = :zos_fts) do c
    masked!(c, time_mean_drop(c, :zos_fts))
end

LOADERS[:heat_flux] = disk_cached(:heat_flux; source_fts_syms = :hfds_fts) do c
    masked!(c, time_mean_drop(c, :hfds_fts) .* (ρ_ocean * cp_ocean))
end
LOADERS[:freshwater_flux] = disk_cached(:freshwater_flux; source_fts_syms = :wfo_fts) do c
    masked!(c, time_mean_drop(c, :wfo_fts))
end

LOADERS[:ssh_squared_mean] = disk_cached(:ssh_squared_mean; source_fts_syms = :zossq_fts) do c
    time_mean_drop(c, :zossq_fts)
end
# `:ssh_variance` is just a subtraction of two cached arrays — keep in memory only.
LOADERS[:ssh_variance] = c -> masked!(c, get_field(c, :ssh_squared_mean) .- get_field(c, :ssh) .^ 2)

# Wind stress is stored as kinematic flux; flip sign (atm→ocean
# downward = positive CMIP convention) and multiply by ρ to get N/m².
# Then bring τx/τy from cell faces down to centers if the grid stores
# them staggered. Cache both components together in a single pair so
# the time-mean of `tauuo`/`tauvo` runs once even if both stresses
# are requested.
centered_stress(c, fts_sym) = -time_mean_drop(c, fts_sym) .* ρ_ocean

LOADERS[:wind_stress_pair] = disk_cached(:wind_stress_pair; source_fts_syms = :tauuo_fts) do c
    τx = centered_stress(c, :tauuo_fts)
    τy = centered_stress(c, :tauvo_fts)
    size(τx, 1) == size(τy, 1) + 1 && (τx = (τx[1:end-1, :] .+ τx[2:end, :]) ./ 2)
    size(τy, 2) == size(τx, 2) + 1 && (τy = (τy[:, 1:end-1] .+ τy[:, 2:end]) ./ 2)
    return (masked!(c, τx), masked!(c, τy))
end

LOADERS[:zonal_wind_stress]      = c -> get_field(c, :wind_stress_pair)[1]
LOADERS[:meridional_wind_stress] = c -> get_field(c, :wind_stress_pair)[2]

#####
##### Near-surface currents (top model cell, rotated to geographic E/N)
#####
#
# The (u, v) components stored by the model are aligned with the grid axes,
# which on ORCA / tripolar grids do NOT coincide with geographic east/north
# above ~60°N. Regridding the grid-aligned components directly to lat-lon
# would mix vectors that point in different physical directions, so we
# rotate to (uE, vN) on the model grid first and then let the existing
# conservative regridder produce the lat-lon maps.

# Rotation between the grid-i direction and geographic east at cell
# centers. Computed from the longitude/latitude of the two F-C-C nodes
# bracketing each cell; `λnode`/`φnode` already handle periodic-x halos,
# so `i = Nx → Nx+1` wraps correctly.
function compute_grid_rotation_centers(grid)
    Nx, Ny, _ = size(grid)
    F, C = Face(), Center()
    cosθ = Array{Float64}(undef, Nx, Ny)
    sinθ = Array{Float64}(undef, Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        λw = λnode(i,   j, 1, grid, F, C, C)
        λe = λnode(i+1, j, 1, grid, F, C, C)
        φw = φnode(i,   j, 1, grid, F, C, C)
        φe = φnode(i+1, j, 1, grid, F, C, C)
        φ  = φnode(i,   j, 1, grid, C, C, C)
        Δλ = mod(λe - λw + 180, 360) - 180
        Δφ = φe - φw
        dE = cosd(φ) * Δλ
        dN = Δφ
        r  = hypot(dE, dN)
        cosθ[i, j] = r > 0 ? dE / r : 1.0
        sinθ[i, j] = r > 0 ? dN / r : 0.0
    end
    return (cosθ, sinθ)
end

LOADERS[:grid_rotation_pair] = c -> compute_grid_rotation_centers(get_field(c, :grid))

# Time-mean u, v at the top model cell (k = Nz), interpolated from their
# native (F, C, C) / (C, F, C) locations to cell centers via `@at` so the
# tripolar fold is handled correctly, then rotated from grid-aligned
# (u, v) to geographic (uE, vN). Cached together so the two FTS reads
# happen once even if both components are requested.
LOADERS[:near_surface_velocity_pair] = disk_cached(:near_surface_velocity_pair; source_fts_syms = :uo_fts) do c
    grid = get_field(c, :grid)
    Nx, Ny, Nz = size(grid)
    u_fts = get_field(c, :uo_fts)
    v_fts = get_field(c, :vo_fts)

    u = Field{Face,   Center, Center}(grid)
    v = Field{Center, Face,   Center}(grid)
    u_cc = Field(@at((Center, Center, Center), u))
    v_cc = Field(@at((Center, Center, Center), v))

    idx = in_window(u_fts; start_time = c.start_time, stop_time = c.stop_time)
    isempty(idx) && error("No snapshots in [$(c.start_time), $(c.stop_time)]")

    u_sum = zeros(Nx, Ny)
    v_sum = zeros(Nx, Ny)
    for n in idx
        set!(u, u_fts[n])
        set!(v, v_fts[n])
        compute!(u_cc)
        compute!(v_cc)
        u_sum .+= Array(view(interior(u_cc), :, :, Nz))
        v_sum .+= Array(view(interior(v_cc), :, :, Nz))
    end
    ū = u_sum ./ length(idx)
    v̄ = v_sum ./ length(idx)

    cosθ, sinθ = get_field(c, :grid_rotation_pair)
    uE = ū .* cosθ .- v̄ .* sinθ
    vN = ū .* sinθ .+ v̄ .* cosθ
    return (masked!(c, uE), masked!(c, vN))
end

LOADERS[:near_surface_zonal_velocity]      = c -> get_field(c, :near_surface_velocity_pair)[1]
LOADERS[:near_surface_meridional_velocity] = c -> get_field(c, :near_surface_velocity_pair)[2]

#####
##### Equatorial-undercurrent (EUC) sections
#####
#
# Two section diagnostics modeled on Fig. 5 of Ringler et al. (2013):
#   * Meridional section of zonal velocity at 140°W, lat ∈ [-8°, 10°]
#   * Zonal section along the equator from 143°E to 95°W
# Both are extracted directly on the model grid (ORCA / tripolar grids are
# essentially lat-lon south of ~60°N, so no horizontal regridding is
# needed in the tropical Pacific). Velocities are rotated to geographic
# east per cell — a no-op in the un-rotated part of the grid, but kept
# for consistency with the surface-current pipeline.
#
# Observations (Johnson et al. 2002) are stubbed for now; the loaders
# `:euc_obs_meridional` / `:euc_obs_equatorial` return `nothing` so the
# figure can render model-only.

# Smallest signed angular distance |λ₁ − λ₂| on the periodic-360 circle.
@inline lon_distance(λ₁, λ₂) = abs(mod(λ₁ - λ₂ + 180, 360) - 180)

function nearest_i(grid, j, target_lon)
    Nx = size(grid, 1)
    best_i, best_d = 1, Inf
    for i in 1:Nx
        λ = mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360)
        d = lon_distance(λ, mod(target_lon, 360))
        d < best_d && (best_d = d; best_i = i)
    end
    return best_i
end

function nearest_j(grid, i, target_lat)
    Ny = size(grid, 2)
    best_j, best_d = 1, Inf
    for j in 1:Ny
        φ = φnode(i, j, 1, grid, Center(), Center(), Center())
        d = abs(φ - target_lat)
        d < best_d && (best_d = d; best_j = j)
    end
    return best_j
end

# Returns (is, js, lats): for each j in the latitude window, the i index
# closest to `target_lon` and the corresponding cell-center latitude.
function meridional_section_path(grid, target_lon; lat_range = (-8.0, 10.0))
    Ny = size(grid, 2)
    is   = Int[]
    js   = Int[]
    lats = Float64[]
    for j in 1:Ny
        i = nearest_i(grid, j, target_lon)
        φ = φnode(i, j, 1, grid, Center(), Center(), Center())
        if lat_range[1] <= φ <= lat_range[2]
            push!(is, i); push!(js, j); push!(lats, φ)
        end
    end
    return is, js, lats
end

# Returns (is, j_eq, lons): a single equatorial j and the i indices whose
# cell-center longitude (in 0..360) falls in `lon_range`.
function equatorial_section_path(grid; lon_range = (143.0, 265.0))
    Nx = size(grid, 1)
    j_eq = nearest_j(grid, nearest_i(grid, size(grid, 2) ÷ 2, 180.0), 0.0)
    is   = Int[]
    lons = Float64[]
    for i in 1:Nx
        λ = mod(λnode(i, j_eq, 1, grid, Center(), Center(), Center()), 360)
        if lon_range[1] <= λ <= lon_range[2]
            push!(is, i); push!(lons, λ)
        end
    end
    return is, j_eq, lons
end

# Streaming extraction: per snapshot, interpolate u (F,C,C) and v (C,F,C)
# to cell centers via `@at`, then accumulate only along the two section
# paths. Avoids materializing a 3-D time-mean (which would be many GB on
# tenth-degree grids). After averaging, rotate to geographic east using
# the 2-D `(cosθ, sinθ)` at each cell. Disk-cached; sections themselves
# are tiny so the cache stays small even at high resolution.
LOADERS[:euc_sections] = disk_cached(:euc_sections; source_fts_syms = :uo_fts) do c
    grid = get_field(c, :grid)
    Nx, Ny, Nz = size(grid)
    cosθ, sinθ = get_field(c, :grid_rotation_pair)
    u_fts = get_field(c, :uo_fts)
    v_fts = get_field(c, :vo_fts)

    is_m,  js_m,  lats_m   = meridional_section_path(grid, -140.0; lat_range = (-8.0, 10.0))
    is_eq, j_eq,  lons_eq  = equatorial_section_path(grid; lon_range = (143.0, 265.0))

    u    = Field{Face,   Center, Center}(grid)
    v    = Field{Center, Face,   Center}(grid)
    u_cc = Field(@at((Center, Center, Center), u))
    v_cc = Field(@at((Center, Center, Center), v))

    Nm  = length(js_m)
    Neq = length(is_eq)
    u_m_sum  = zeros(Nm,  Nz)
    v_m_sum  = zeros(Nm,  Nz)
    u_eq_sum = zeros(Neq, Nz)
    v_eq_sum = zeros(Neq, Nz)

    idx = in_window(u_fts; start_time = c.start_time, stop_time = c.stop_time)
    isempty(idx) && error("No snapshots in [$(c.start_time), $(c.stop_time)]")
    for n in idx
        set!(u, u_fts[n]); compute!(u_cc)
        set!(v, v_fts[n]); compute!(v_cc)
        ui = interior(u_cc)
        vi = interior(v_cc)
        @inbounds for jdx in 1:Nm
            i, j = is_m[jdx], js_m[jdx]
            for k in 1:Nz
                u_m_sum[jdx, k] += ui[i, j, k]
                v_m_sum[jdx, k] += vi[i, j, k]
            end
        end
        @inbounds for ieq in 1:Neq
            i = is_eq[ieq]
            for k in 1:Nz
                u_eq_sum[ieq, k] += ui[i, j_eq, k]
                v_eq_sum[ieq, k] += vi[i, j_eq, k]
            end
        end
    end

    Nt = length(idx)
    uE_m  = similar(u_m_sum)
    uE_eq = similar(u_eq_sum)
    @inbounds for jdx in 1:Nm
        i, j = is_m[jdx], js_m[jdx]
        cθ = cosθ[i, j]; sθ = sinθ[i, j]
        for k in 1:Nz
            uE_m[jdx, k] = (u_m_sum[jdx, k] * cθ - v_m_sum[jdx, k] * sθ) / Nt
        end
    end
    @inbounds for ieq in 1:Neq
        i = is_eq[ieq]
        cθ = cosθ[i, j_eq]; sθ = sinθ[i, j_eq]
        for k in 1:Nz
            uE_eq[ieq, k] = (u_eq_sum[ieq, k] * cθ - v_eq_sum[ieq, k] * sθ) / Nt
        end
    end

    depth = collect(znodes(grid, Center()))
    return (uE_meridional = uE_m,    lats_meridional = lats_m,
            uE_equatorial = uE_eq,   lons_equatorial = lons_eq,
            depth = depth)
end

LOADERS[:euc_meridional_section] = c -> get_field(c, :euc_sections).uE_meridional
LOADERS[:euc_meridional_lats]    = c -> get_field(c, :euc_sections).lats_meridional
LOADERS[:euc_equatorial_section] = c -> get_field(c, :euc_sections).uE_equatorial
LOADERS[:euc_equatorial_lons]    = c -> get_field(c, :euc_sections).lons_equatorial
LOADERS[:euc_depth]              = c -> get_field(c, :euc_sections).depth

# Johnson et al. (2002) observational climatology — stub. Returns `nothing`
# so the figure can render model-only; replace with a real loader when
# the obs file is available locally or via a verified download URL.
LOADERS[:euc_obs_meridional] = c -> nothing
LOADERS[:euc_obs_equatorial] = c -> nothing

#####
##### Sea-ice concentration (mean + monthly bins)
#####

LOADERS[:sic_mean_and_monthly] = disk_cached(:sic_mean_and_monthly; source_fts_syms = :sic_fts) do c
    compute_mean_and_monthly(get_field(c, :sic_fts); start_time = c.start_time, stop_time = c.stop_time)
end

function sic_month(c, m)
    _, monthly = get_field(c, :sic_mean_and_monthly)
    isnothing(monthly[m]) ? nothing : masked!(c, dropdims(monthly[m]; dims = 3))
end

LOADERS[:sic_mean]      = c -> masked!(c, dropdims(get_field(c, :sic_mean_and_monthly)[1]; dims = 3))
LOADERS[:sic_march]     = c -> sic_month(c, 3)
LOADERS[:sic_september] = c -> sic_month(c, 9)

# Observed monthly SIC climatology (HadISST1, 1979–2007) lives on the
# shared 1° lat-lon grid by construction (see `hadisst_sic_climatology`
# in common.jl). It is case-independent, but exposed through each
# per-case cache so the same `get_field(c, :sic_*_bias_latlon)` pattern
# composes uniformly with the other lat-lon loaders.
LOADERS[:hadisst_sic_monthly] = c -> hadisst_sic_climatology()

function hadisst_month_latlon(c, m)
    monthly = get_field(c, :hadisst_sic_monthly)
    return isnothing(monthly) ? nothing : monthly[m]
end

LOADERS[:sic_march_obs_latlon]     = c -> hadisst_month_latlon(c, 3)
LOADERS[:sic_september_obs_latlon] = c -> hadisst_month_latlon(c, 9)

# Model minus obs SIC bias on the shared 1° lat-lon grid. `nothing` in
# from either side → `nothing` out so fig06 degrades gracefully if the
# HadISST download fails or the case lacks the corresponding monthly bin.
#
# NaN propagates: if either model or obs is NaN at a cell, the bias is
# NaN there. This is intentional — at 1° regridded resolution the
# model's land mask and HadISST's land mask disagree on a thin strip of
# coastal cells, and forcing those cells to 0 would produce spurious
# `model − 0 = model` or `0 − obs = −obs` values that contaminate the
# bias panels with large near-coastal artefacts. Letting them be NaN
# renders them as no-data (the Natural-Earth land polygon underneath
# the contourf hides them where it can; the residual coastal-strip
# pixels stay grey, correctly signalling "no comparison possible").
function sic_bias_latlon(c, model_sym, obs_sym)
    model = get_field(c, model_sym)
    obs   = get_field(c, obs_sym)
    (isnothing(model) || isnothing(obs)) && return nothing
    return model .- obs
end

LOADERS[:sic_march_bias_latlon]     = c -> sic_bias_latlon(c, :sic_march_latlon,     :sic_march_obs_latlon)
LOADERS[:sic_september_bias_latlon] = c -> sic_bias_latlon(c, :sic_september_latlon, :sic_september_obs_latlon)

#####
##### Mixed-layer depth (kernel diagnostic, monthly min/max)
#####

LOADERS[:mld_monthly] = disk_cached(:mld_monthly; source_fts_syms = :mld_fts) do c
    compute_monthly_means(get_field(c, :mld_fts); start_time = c.start_time, stop_time = c.stop_time)
end

# Stream a per-cell reduction across all available monthly bins of
# `sym`, avoiding the ~Nx·Ny·12 temp cube `cat(...; dims = 3)` would
# allocate. Initializes from the first available month, then folds the
# rest in place via `reduce.(acc, monthly[m])`.
function reduce_monthly(c, sym, reduce)
    monthly = get_field(c, sym)
    avail   = findall(!isnothing, monthly)
    isempty(avail) && error("$sym has no available monthly bins")
    acc = copy(dropdims(monthly[first(avail)]; dims = 3))
    for m in @view avail[2:end]
        acc .= reduce.(acc, dropdims(monthly[m]; dims = 3))
    end
    return acc
end

LOADERS[:mld_min] = c -> masked!(c, reduce_monthly(c, :mld_monthly, min))
LOADERS[:mld_max] = c -> masked!(c, reduce_monthly(c, :mld_monthly, max))

#####
##### WOA temperature & salinity on case grid (TEOS-10)
#####
#
# WOA `t_an` is in-situ temperature (°C) and `s_an` is Practical Salinity
# (PSS-78), but the OMIP ocean prognostic state is TEOS-10 Conservative
# Temperature (Θ) and Absolute Salinity (S_A). Comparing the model
# directly against the raw WOA fields produces a near-constant
# +0.165 g/kg salinity bias (the S_A − S_P offset) and a smaller
# salinity-dependent temperature bias. The pair loader interpolates
# both fields onto the case grid and runs `woa_to_teos10!` once so
# every downstream loader (`:sst_bias`, `:sss_bias`, `:zonal_woa_*`,
# `:zonal_*_bias`) compares apples to apples.

function woa_teos10_pair(c)
    grid = get_field(c, :grid)
    woaT = Field(Metadatum(:temperature; dataset = WOAAnnual()), CPU())
    woaS = Field(Metadatum(:salinity;    dataset = WOAAnnual()), CPU())
    T = CenterField(grid)
    S = CenterField(grid)
    interpolate!(T, woaT)
    interpolate!(S, woaS)
    return (Array(interior(T)), Array(interior(S)))
end

LOADERS[:woa_pair]        = woa_teos10_pair
LOADERS[:woa_temperature] = c -> get_field(c, :woa_pair)[1]
LOADERS[:woa_salinity]    = c -> get_field(c, :woa_pair)[2]

LOADERS[:sst_bias] = c -> masked!(c, get_field(c, :sst) .- get_field(c, :woa_temperature)[:, :, end])
LOADERS[:sss_bias] = c -> masked!(c, get_field(c, :sss) .- get_field(c, :woa_salinity)[:, :, end])

#####
##### ECCO SSH
#####

LOADERS[:ecco_ssh] = c -> masked!(c, ecco_ssh_on_grid(get_field(c, :grid)))

LOADERS[:ssh_bias_ecco] = c -> masked!(c, begin
    η     = get_field(c, :ssh)
    η_ref = get_field(c, :ecco_ssh)
    (η .- mean(filter(isfinite, η))) .- (η_ref .- mean(filter(isfinite, η_ref)))
end)

#####
##### dBM mixed-layer climatology
#####

LOADERS[:dbm_mld_monthly] = c -> dbm_mld_climatology_on_grid(get_field(c, :grid))

dbm_extreme(c, reduce) = masked!(c, let m = get_field(c, :dbm_mld_monthly)
    isnothing(m) ? nothing : dropdims(reduce(m; dims = 3); dims = 3)
end)

LOADERS[:mld_min_dbm] = c -> dbm_extreme(c, minimum)
LOADERS[:mld_max_dbm] = c -> dbm_extreme(c, maximum)

#####
##### NCEP wind-stress climatology
#####

LOADERS[:ncep_wind_stress_pair]            = c -> ncep_wind_stress_on_grid(get_field(c, :grid))
LOADERS[:ncep_zonal_wind_stress]           = c -> masked!(c, get_field(c, :ncep_wind_stress_pair)[1])
LOADERS[:ncep_meridional_wind_stress]      = c -> masked!(c, get_field(c, :ncep_wind_stress_pair)[2])
LOADERS[:zonal_wind_stress_bias_ncep]      = c -> maybe_diff!(c, :zonal_wind_stress,      :ncep_zonal_wind_stress)
LOADERS[:meridional_wind_stress_bias_ncep] = c -> maybe_diff!(c, :meridional_wind_stress, :ncep_meridional_wind_stress)

#####
##### Sea-ice integrals (heavy: per-snapshot loop over sithick + sic)
#####

LOADERS[:sea_ice_diagnostics] = disk_cached(:sea_ice_diagnostics; source_fts_syms = :sithick_fts) do c
    compute_ice_diagnostics(c.run_dir, c.prefix, get_field(c, :grid);
                            start_time = c.start_time,
                            stop_time  = c.stop_time)
end

#####
##### Time-series scalars + horizontal-mean profiles
#####

LOADERS[:global_mean_temperature_timeseries] =
    disk_cached(:global_mean_temperature_timeseries; source_fts_syms = :tosga_fts) do c
    total_jld2_scalar_timeseries(get_field(c, :averages_file), "tosga")
end
LOADERS[:global_mean_salinity_timeseries] =
    disk_cached(:global_mean_salinity_timeseries; source_fts_syms = :soga_fts) do c
    total_jld2_scalar_timeseries(get_field(c, :averages_file), "soga")
end

LOADERS[:time_in_years] = c ->
    total_jld2_timeseries_times(get_field(c, :averages_file)) ./ (365.25 * 24 * 3600)

LOADERS[:horizontal_mean_temperature_profile] = c ->
    vec(compute_time_mean(get_field(c, :to_h_fts);
                           start_time = c.start_time,
                           stop_time  = c.stop_time))

LOADERS[:horizontal_mean_salinity_profile] = c ->
    vec(compute_time_mean(get_field(c, :so_h_fts);
                           start_time = c.start_time,
                           stop_time  = c.stop_time))

# Per-snapshot horizontal-mean drift relative to the first snapshot.
# Intentionally ignores `c.start_time`/`c.stop_time` — drift is always
# referenced to t = 0 of the run, not to the case averaging window
# (figures 15 and 20 plot the full record so the spinup is visible).
function profile_drift(c, fts_sym)
    fts = get_field(c, fts_sym)
    Nt  = length(fts.times)
    Nz  = size(fts[1], 3)
    Δ   = zeros(Nt, Nz)
    for n in 1:Nt
        Δ[n, :] .= vec(interior(fts[n]))
    end
    Δ .-= reshape(Δ[1, :], 1, :)
    return Δ
end

LOADERS[:temperature_drift] = disk_cached(:temperature_drift; source_fts_syms = :to_h_fts) do c
    profile_drift(c, :to_h_fts)
end
LOADERS[:salinity_drift] = disk_cached(:salinity_drift; source_fts_syms = :so_h_fts) do c
    profile_drift(c, :so_h_fts)
end

#####
##### Global-mean kinetic energy from u, v snapshots
#####

# Per-cell volume-weighted kinetic energy at cell centers, written as a
# `KernelFunctionOperation` so `compute!` runs as a single fused kernel
# — much faster than the previous
# `Field(@at((Center, Center, Center), u*u + v*v))` AbstractOperation
# tree, which built intermediate buffers for `u*u`, `v*v`, and the
# interpolations to centers.
#
# `ψ²` squares its argument before the surrounding `ℑxᶜᵃᵃ` / `ℑyᵃᶜᵃ`
# average to centers, matching the original "square then interpolate"
# semantics of the `@at` expression. (Same pattern as `ϕ²` /
# `τᶜᶜᶜ` in `friction_velocity.jl`.) Multiplying by `Vᶜᶜᶜ` here gives
# the cell's KE contribution, so `sum(field) / sum(V) ` is the
# volume-averaged KE — and on an `ImmersedBoundaryGrid` both sums
# automatically skip immersed cells, replacing the old hand-rolled
# `ocean_mask_3d` book-keeping.
@inline ψ²(i, j, k, grid, ψ) = @inbounds ψ[i, j, k]^2

@inline uu_plus_vv_ccc(i, j, k, grid, u, v) =
    Vᶜᶜᶜ(i, j, k, grid) *
    (ℑxᶜᵃᵃ(i, j, k, grid, ψ², u) + ℑyᵃᶜᵃ(i, j, k, grid, ψ², v)) / 2


LOADERS[:kinetic_energy_pair] = disk_cached(:kinetic_energy_pair; source_fts_syms = :uo_fts) do c
    grid   = get_field(c, :grid)
    u_fts  = get_field(c, :uo_fts)
    v_fts  = get_field(c, :vo_fts)

    u  = Field{Face,   Center, Center}(grid)
    v  = Field{Center, Face,   Center}(grid)
    e_op = KernelFunctionOperation{Center, Center, Center}(uu_plus_vv_ccc, grid, u, v)
    V_op = KernelFunctionOperation{Center, Center, Center}(Vᶜᶜᶜ, grid)
    e    = Field(e_op)
    V    = sum(V_op)

    Nt = length(u_fts.times)
    ke = zeros(Nt)
    for n in 1:Nt
        set!(u, u_fts[n])
        set!(v, v_fts[n])
        compute!(e)
        ke[n] = sum(e) / V
    end
    return (ke, u_fts.times ./ (365.25 * 24 * 3600))
end

LOADERS[:kinetic_energy]               = c -> get_field(c, :kinetic_energy_pair)[1]
LOADERS[:kinetic_energy_time_in_years] = c -> get_field(c, :kinetic_energy_pair)[2]

#####
##### Zonal means
#####
#
# Target lat-lon grid is shared across all cases.

const ZONAL_NLON = 360
const ZONAL_NLAT = 180

# Target grid + destination field are tiny and shared across cases, so
# bind them eagerly at module load. Type-stable; consumers (the zonal
# regridder loader) are themselves lazy via the LOADERS registry.
const ZONAL_LATLON_GRID = LatitudeLongitudeGrid(CPU();
    size      = (ZONAL_NLON, ZONAL_NLAT, 1),
    longitude = (0, 360), latitude = (-90, 90), z = (0, 1))

const ZONAL_LATLON_DST_FIELD = Field{Center, Center, Nothing}(ZONAL_LATLON_GRID)

# Cell centers on the shared lat-lon grid. Used as axis coordinates by
# `surface_panel!` so every regridded map renders with physical lon/lat
# axes regardless of the underlying model grid (ORCA, tripolar, …).
const LATLON_LON_CENTERS = collect(λnodes(ZONAL_LATLON_GRID, Center()))
const LATLON_LAT_CENTERS = collect(φnodes(ZONAL_LATLON_GRID, Center()))

zonal_latitude_centers() = collect(φnodes(ZONAL_LATLON_GRID, Center()))

# Conservatively regrid each level of `data_3d` (weighted by
# `ocean_mask_3d`) onto the lat-lon target grid, then average per latitude row.
function compute_zonal_mean(data_3d, ocean_mask_3d, regridder, Nlon, Nlat)
    Nz    = size(data_3d, 3)
    zonal = fill(NaN, Nlat, Nz)
    fdata = zeros(Nlon * Nlat)
    fmask = zeros(Nlon * Nlat)
    A     = regridder.dst_areas
    for k in 1:Nz
        ConservativeRegridding.regrid!(fdata, regridder, vec(data_3d[:, :, k] .* ocean_mask_3d[:, :, k]))
        ConservativeRegridding.regrid!(fmask, regridder, vec(ocean_mask_3d[:, :, k]))
        Wd = reshape(fdata .* A, Nlon, Nlat)
        Wm = reshape(fmask .* A, Nlon, Nlat)
        for j in 1:Nlat
            m = sum(@view Wm[:, j])
            m > 0 && (zonal[j, k] = sum(@view Wd[:, j]) / m)
        end
    end
    return zonal
end

# Session-level cache of `ConservativeRegridding.Regridder` keyed by a
# content fingerprint of the source grid. Two cases that share the same
# physical grid (typical: many cases run against the same config, e.g.
# ORCA1) build the regridder exactly once and then share it. The
# `Regridder` constructor takes minutes per case at 1/10° because it
# allocates sparse weight matrices over O(N²) candidate cell pairs, so
# this cache saves wall-clock proportional to the case count.
const REGRIDDER_CACHE = Dict{Any, Any}()

# Fingerprint identifying "the same grid for regridding purposes":
# type + size + bottom-topography hash. Bottom topography captures the
# entire land mask (which is what `ConservativeRegridding` actually
# depends on beyond the lat/lon axes), and for a fixed config grid the
# axes are uniquely determined by type+size. If two grids agree on all
# three they produce bit-identical regridders.
function regridder_cache_key(grid)
    key = Any[string(nameof(typeof(grid))), size(grid)]
    if grid isa ImmersedBoundaryGrid
        bh = Array(interior(grid.immersed_boundary.bottom_height, :, :, 1))
        push!(key, hash(bh))
    end
    return Tuple(key)
end

LOADERS[:zonal_regridder] = c -> begin
    grid = get_field(c, :grid)
    key  = regridder_cache_key(grid)
    cached = get(REGRIDDER_CACHE, key, nothing)
    if !isnothing(cached)
        @info "  Reusing zonal regridder for $(c.label) (cache hit on grid fingerprint)"
        return cached
    end
    src = Field{Center, Center, Nothing}(grid)
    @info "  Building zonal regridder for $(c.label) (may take a few minutes)..."
    rg = ConservativeRegridding.Regridder(ZONAL_LATLON_DST_FIELD, src; progress = true)
    REGRIDDER_CACHE[key] = rg
    return rg
end

# Surface regridder is the *same* (model 2-D grid → 1° lat-lon) regridder
# used for zonal means. Alias rather than rebuild so each case pays the
# regridder-construction cost exactly once.
LOADERS[:surface_regridder] = c -> get_field(c, :zonal_regridder)

"""
    regrid_surface_to_latlon(c, data_2d)

Conservatively regrid a 2-D surface field `data_2d` (defined on the case's
native ocean grid) onto the shared 1° lat-lon grid. Returns a
`(NLON, NLAT)` array, NaN-masked where the destination cell has zero
ocean coverage. `nothing` in → `nothing` out so it composes through
optional climatology fields.
"""
function regrid_surface_to_latlon(c::CaseCache, data_2d)
    isnothing(data_2d) && return nothing
    regridder    = get_field(c, :surface_regridder)
    surface_mask = Float64.(get_field(c, :ocean_mask_3d)[:, :, end])

    fdata = zeros(ZONAL_NLON * ZONAL_NLAT)
    fcov  = zeros(ZONAL_NLON * ZONAL_NLAT)
    clean = replace(data_2d, NaN => zero(eltype(data_2d)))

    ConservativeRegridding.regrid!(fdata, regridder, vec(clean .* surface_mask))
    ConservativeRegridding.regrid!(fcov,  regridder, vec(surface_mask))

    raw = reshape(fdata, ZONAL_NLON, ZONAL_NLAT)
    cov = reshape(fcov,  ZONAL_NLON, ZONAL_NLAT)
    out = fill(NaN, ZONAL_NLON, ZONAL_NLAT)
    @inbounds for i in eachindex(cov)
        cov[i] > 0 && (out[i] = raw[i] / cov[i])
    end
    return out
end

# ssh_variance is η² mean minus η_mean²; SSH RMS = √(variance), clipped
# at zero to guard against tiny negative values from finite-sample noise.
LOADERS[:ssh_rms] = c -> begin
    v = get_field(c, :ssh_variance)
    isnothing(v) && return nothing
    return sqrt.(max.(v, zero(eltype(v))))
end

# Build a `:foo_latlon` loader for every surface field used by figures.
# Each is a lazy lat-lon regrid of the existing case-grid loader, cached
# in-memory only (no disk persistence — the case-grid version already has
# disk caching where appropriate, and regridding is fast compared to
# rebuilding the regridder).
const SURFACE_LATLON_FIELDS = (
    :sst, :sss, :ssh,
    :sst_bias, :sss_bias, :ssh_bias_ecco,
    :heat_flux, :freshwater_flux,
    :ssh_rms,
    :zonal_wind_stress, :meridional_wind_stress,
    :zonal_wind_stress_bias_ncep, :meridional_wind_stress_bias_ncep,
    :near_surface_zonal_velocity, :near_surface_meridional_velocity,
    :sic_mean, :sic_march, :sic_september,
    :mld_min, :mld_max, :mld_min_dbm, :mld_max_dbm,
)

for sym in SURFACE_LATLON_FIELDS
    let s = sym
        LOADERS[Symbol(s, :_latlon)] = c -> regrid_surface_to_latlon(c, get_field(c, s))
    end
end

# Near-surface current speed on the lat-lon grid: take the magnitude of
# the regridded geographic components rather than regridding |u| from the
# model grid. This gives the speed of the (area-mean) Eulerian current
# instead of the area-mean of the speed, which is the conventional
# "mean current" quantity (and consistent with the model fields it's
# compared to).
LOADERS[:near_surface_speed_latlon] = c -> hypot.(get_field(c, :near_surface_zonal_velocity_latlon),
                                                  get_field(c, :near_surface_meridional_velocity_latlon))

"""
    surface_panel!(fig, pos, data; projection = "+proj=robin", kwargs...)

Plot a regridded surface field as a filled-contour map in `projection`
(Robinson by default). Coastlines and land polygons (Natural Earth)
are drawn automatically.

`data` must be `(NLON, NLAT)` on the shared 1° lat-lon grid — i.e.
produced by the matching `:foo_latlon` loader or
`regrid_surface_to_latlon`. All `geo_panel!` kwargs are forwarded.
"""
function surface_panel!(fig, pos, data; projection = "+proj=robin", kwargs...)
    return geo_panel!(fig, pos, data;
                      x = LATLON_LON_CENTERS,
                      y = LATLON_LAT_CENTERS,
                      projection,
                      kwargs...)
end

"""
    polar_panel!(fig, pos, data; hemisphere = :north, lat_cutoff = 45.0, kwargs...)

Plot a regridded surface field with a polar stereographic projection,
clipped to the requested hemisphere (`:north` or `:south`) and to
latitudes poleward of `lat_cutoff` (degrees). Used by the sea-ice
concentration figure to show NH and SH separately.

`data` must be `(NLON, NLAT)` on the shared 1° lat-lon grid; the
function extracts only the latitude band of interest and computes
projected limits (via Proj) so the polar cap fills the panel instead
of being squashed into the centre of a full-globe stereographic view.
"""
function polar_panel!(fig, pos, data;
                      hemisphere = :north, lat_cutoff = 45.0,
                      obs_contour = nothing, kwargs...)
    proj = hemisphere == :north ?
        "+proj=stere +lat_0=90 +lat_ts=70" :
        "+proj=stere +lat_0=-90 +lat_ts=-70"
    j_keep = hemisphere == :north ?
        findall(lat -> lat >= lat_cutoff,  LATLON_LAT_CENTERS) :
        findall(lat -> lat <= -lat_cutoff, LATLON_LAT_CENTERS)
    lat_sub     = LATLON_LAT_CENTERS[j_keep]
    data_sub    = data[:, j_keep]
    obs_sub     = isnothing(obs_contour) ? nothing : obs_contour[:, j_keep]
    latlims     = hemisphere == :north ? (lat_cutoff, 90.0) : (-90.0, -lat_cutoff)
    return geo_panel!(fig, pos, data_sub;
                      x = LATLON_LON_CENTERS,
                      y = lat_sub,
                      projection = proj,
                      latlims = latlims,
                      lonlims = (-180.0, 180.0),
                      obs_contour = obs_sub,
                      kwargs...)
end

# 3-D time-means (loaded on demand for zonal-section figures). Not
# disk-cached because each is several hundred megabytes per ORCA case;
# the smaller `:zonal_*` derivatives below are.
for (sym, fts) in ((:time_mean_temperature_3d, :to_fts),
                   (:time_mean_salinity_3d,    :so_fts),
                   (:time_mean_buoyancy_3d,    :bo_fts))
    LOADERS[sym] = let f = fts
        c -> compute_time_mean(get_field(c, f); start_time = c.start_time, stop_time = c.stop_time)
    end
end
LOADERS[:initial_buoyancy_3d] = c -> Array(interior(get_field(c, :bo_fts)[1]))

# Zonal mean of any 3-D field already in the cache.
zonal_of(c, sym) = compute_zonal_mean(get_field(c, sym), get_field(c, :ocean_mask_3d),
                                       get_field(c, :zonal_regridder), ZONAL_NLON, ZONAL_NLAT)

# Each zonal section: (sym, source-3D, source-FTS for invalidation).
# WOA entries pass `nothing` since they don't depend on model output.
for (sym, src, fts) in ((:zonal_temperature,       :time_mean_temperature_3d, :to_fts),
                        (:zonal_salinity,          :time_mean_salinity_3d,    :so_fts),
                        (:zonal_buoyancy,          :time_mean_buoyancy_3d,    :bo_fts),
                        (:zonal_initial_buoyancy,  :initial_buoyancy_3d,      :bo_fts),
                        (:zonal_woa_temperature,   :woa_temperature,          nothing),
                        (:zonal_woa_salinity,      :woa_salinity,             nothing))
    LOADERS[sym] = let s = src
        wrapped = c -> zonal_of(c, s)
        isnothing(fts) ? disk_cached(wrapped, sym) :
                         disk_cached(wrapped, sym; source_fts_syms = fts)
    end
end

LOADERS[:zonal_temperature_bias] = c -> get_field(c, :zonal_temperature) .-
                                          get_field(c, :zonal_woa_temperature)
LOADERS[:zonal_salinity_bias]    = c -> get_field(c, :zonal_salinity) .-
                                          get_field(c, :zonal_woa_salinity)
LOADERS[:zonal_buoyancy_drift]   = c -> get_field(c, :zonal_buoyancy) .-
                                          get_field(c, :zonal_initial_buoyancy)

# Zonal MLD: regrid the 2-D surface MLD field, weighted by the surface
# ocean mask. NaNs (land) become 0 so they don't poison the regrid.
function zonal_mld(c, sym)
    raw = get_field(c, sym)
    isnothing(raw) && return nothing
    raw_3d  = reshape(replace(raw, NaN => zero(eltype(raw))), size(raw)..., 1)
    surface = get_field(c, :ocean_mask_3d)[:, :, end:end]
    return vec(compute_zonal_mean(raw_3d, surface, get_field(c, :zonal_regridder),
                                   ZONAL_NLON, ZONAL_NLAT))
end

LOADERS[:zonal_mld_min]     = c -> zonal_mld(c, :mld_min)
LOADERS[:zonal_mld_max]     = c -> zonal_mld(c, :mld_max)
LOADERS[:zonal_mld_min_dbm] = c -> zonal_mld(c, :mld_min_dbm)
LOADERS[:zonal_mld_max_dbm] = c -> zonal_mld(c, :mld_max_dbm)

#####
##### AMOC streamfunction (Atlantic basin, no regridding)
#####
#
# ORCA / tripolar grids are essentially lat-lon below ~60°N, so the AMOC
# streamfunction is computed per j-row directly on the model grid:
#
#     ψ_atl(j, k) = ∑_{k' ≤ k}  ∑_{i ∈ Atlantic}  vvol(i, j, k')
#
# `vvol = v · Aʸ` is saved by the model on (Center, Face, Center), so the
# zonal sum already carries the proper Δx and the cumulative-z sum carries
# Δz — no extra grid metrics needed offline. The Atlantic basin mask is
# from `Bathymetry.atlantic_ocean_basin`, computed via flood-fill with
# Cape Agulhas and Drake Passage barriers and a 65°N northern cap.

LOADERS[:atlantic_mask_2d] = disk_cached(:atlantic_mask_2d) do c
    basin = atlantic_ocean_basin(get_field(c, :grid))
    return Bool.(dropdims(Array(interior(basin.mask)); dims = 3))
end

# Latitude per j-row at the v-face, averaged across i so the curve is
# scalar in j even on tripolar/ORCA. Below ~60°N this is essentially a
# constant latitude row; the average is just a robust projection above.
# `vvol` is at Face-y, which on `RightFaceFolded` tripolar/ORCA has Ny+1
# interior face cells (the topmost face being the fold line), so we
# produce Ny+1 latitudes to match.
LOADERS[:amoc_latitudes] = c -> begin
    grid = get_field(c, :grid)
    Nx, Ny, _ = size(grid)
    return [mean(φnode(i, j, 1, grid, Center(), Face(), Center()) for i in 1:Nx)
            for j in 1:Ny+1]
end

LOADERS[:amoc] = disk_cached(:amoc; source_fts_syms = :vvol_fts) do c
    vvol_mean = compute_time_mean(get_field(c, :vvol_fts);
                                   start_time = c.start_time,
                                   stop_time  = c.stop_time)
    atl = get_field(c, :atlantic_mask_2d)
    Nx, Ny, Nz = size(vvol_mean)
    transport_per_layer = zeros(Ny, Nz)
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        atl[i, j] || continue
        transport_per_layer[j, k] += vvol_mean[i, j, k]
    end
    # ψ(j, z) = -∫_{-H}^{z} ∑_atl v dx dz'.
    # Oceananigans has k=1 at the bottom, so a cumulative sum along k
    # integrates from -H upward; the leading minus flips the sign so a
    # positive value renders the NADW cell as red on a balance colormap.
    return -cumsum(transport_per_layer; dims = 2) ./ 1e6   # Sv
end

#####
##### AMOC at 26.5°N (RAPID-MOCHA latitude)
#####
#
# `:amoc_latitudes` are at v-faces (length Ny+1); the (Ny, Nz) `:amoc`
# array is indexed by cell centers, so we average adjacent face
# latitudes to find the j closest to 26.5°N.

const RAPID_LATITUDE = 26.5

LOADERS[:amoc_26n_j] = c -> begin
    lats = get_field(c, :amoc_latitudes)
    centers = @views (lats[1:end-1] .+ lats[2:end]) ./ 2
    return argmin(abs.(centers .- RAPID_LATITUDE))
end

LOADERS[:amoc_26n_profile] = c -> begin
    ψ = get_field(c, :amoc)               # (Ny, Nz) Sv
    j = get_field(c, :amoc_26n_j)
    return Vector(ψ[j, :])
end

# Per-snapshot ψ_max at 26.5°N: loop vvol over time, sum across the
# Atlantic at the single j-row, cumulatively integrate in z, take max.
# Cheap per snapshot (O(Nx·Nz)) so this covers the full record even on
# tenth-degree grids. Disk-cached on vvol_fts so reruns are instant.
LOADERS[:amoc_max_timeseries] = disk_cached(:amoc_max_timeseries; source_fts_syms = :vvol_fts) do c
    vvol_fts = get_field(c, :vvol_fts)
    atl      = get_field(c, :atlantic_mask_2d)
    j        = get_field(c, :amoc_26n_j)
    Nx, Ny, Nz = size(get_field(c, :grid))
    Nt = length(vvol_fts.times)
    ψ_max = zeros(Nt)
    for n in 1:Nt
        slice = Array(interior(vvol_fts[n]))
        col = zeros(Nz)
        @inbounds for k in 1:Nz, i in 1:Nx
            atl[i, j] || continue
            col[k] += slice[i, j, k]
        end
        ψ_z = -cumsum(col) ./ 1e6   # Sv, sign matches :amoc
        ψ_max[n] = maximum(ψ_z)
    end
    # Convert to decimal calendar year using the JRA55-do epoch the rest
    # of the visualize pipeline assumes (`compute_monthly_means` uses
    # `DateTime(1958, 1, 1)` as the reference).
    year_start = 1958.0
    times_year = year_start .+ vvol_fts.times ./ years
    return (year = times_year, psi_max = ψ_max)
end

#####
##### Strait transports (offline; depends on per-case grid configuration)
#####

function strait_config_for(case)
    haskey(case, :config) && return case.config
    p = lowercase(case.prefix)
    for cfg in (:tenthdegree, :halfdegree, :orca)
        occursin(string(cfg), p) && return cfg
    end
    return nothing
end

LOADERS[:strait_transports] = c -> begin
    config = strait_config_for(c.case)
    if isnothing(config)
        @warn "Cannot infer strait config for case '$(c.label)' — skipping."
        return nothing
    end
    @info "  $(c.label): computing strait transports ($config)..."
    return strait_transports(config, get_field(c, :fields_file);
                              start_time = 0,
                              stop_time  = Inf)
end
