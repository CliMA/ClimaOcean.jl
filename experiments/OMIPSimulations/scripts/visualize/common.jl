# visualize/common.jl
#
# Imports, monkey-patches, FTS_BACKEND, plotting helpers, time-averaging
# utilities, land/ocean masks, climatology readers (ECCO, dBM, NCEP).
#
# Expects the orchestrator to have already defined:
#   - `cases`        :: Vector{NamedTuple}   (the case list)
#   - `output_dir`   :: String                (where PNGs go)
#
# Defines `obs_cache_dir`, `FTS_BACKEND`, `case_colors`, `savefig`, etc.

# ══════════════════════════════════════════════════════════════
# Constants
# ══════════════════════════════════════════════════════════════

const years    = 365 * 24 * 3600
const ρ_ocean  = 1026.0
const cp_ocean = 3991.86795711963

# ══════════════════════════════════════════════════════════════
# Imports
# ══════════════════════════════════════════════════════════════

using CairoMakie
using GeoMakie
import GeoMakie.GeometryBasics
using Statistics
using Dates
using Downloads
using DelimitedFiles
using JLD2
using NCDatasets
using WorldOceanAtlasTools
using Oceananigans
using Oceananigans.Grids: znodes, λnodes, φnodes, λnode, φnode
using Oceananigans.Fields: interpolate!
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, Vᶜᶜᶜ
using ConservativeRegridding
using NumericalEarth
using NumericalEarth.DataWrangling: Metadatum
using NumericalEarth.DataWrangling.WOA: WOAAnnual
using NumericalEarth: ECCO4Monthly
using OMIPSimulations: strait_transports, woa_to_teos10!

# The Oceananigans monkey-patches that used to live here (split-file
# `FieldTimeSeries` support and the matching `set!` extension), plus their
# JLD2 helpers, now live in `src/oceananigans_patches.jl` of the
# OMIPSimulations module and load when `using OMIPSimulations` runs above.
# We re-bind the helper names below so the rest of this script and `cache.jl`
# keep working unchanged.
using OMIPSimulations: jld2_output_part_paths, with_jld2, jld2_parts,
                       memoize_jld2_part,
                       total_jld2_timeseries_snapshot_count,
                       total_jld2_timeseries_times,
                       last_jld2_timeseries_time,
                       total_jld2_serialized_grid,
                       total_jld2_scalar_timeseries,
                       detect_split_file_path,
                       location_types,
                       rebuild_backend_with_path,
                       rebuild_fts_with_path

# ══════════════════════════════════════════════════════════════
# Output directory + obs cache
# ══════════════════════════════════════════════════════════════

mkpath(output_dir)
@info "Figures will be saved to: $output_dir"

const obs_cache_dir = joinpath(output_dir, "obs_cache")
mkpath(obs_cache_dir)

# Shared backend template — `deepcopy(FTS_BACKEND)` for every FieldTimeSeries
# so each one gets its own independent buffer state. `prefetch = false`
# because multiple FTS share the same JLD2 file and `Prefetched` assumes
# sole-reader access.
const FTS_BACKEND = InMemory(10; prefetch = false)

savefig(fig, name) = save(joinpath(output_dir, name), fig)

# ── Presentation themes ──────────────────────────────────────────────
#
# `:light` is Makie's default (white background, black axis chrome).
# `:dark` returns a transparent figure background with white axis chrome,
# colorbar/legend, and titles — intended for layering the saved PNGs onto
# a dark slide. Apply with
#
#     with_theme(presentation_theme(:dark)) do
#         fig04(caches, labels, cases)
#     end
#
# which is what `render_figures(...; theme = :dark)` does for every
# figure in a single sweep.
#
# Note: `geo_panel!` hard-codes `coast_color = :black` / `land_color =
# :lightgray`; those are plot content, not axis chrome, and are not
# touched by the theme. Override at the call site if a given map needs
# them flipped.

function presentation_theme(name::Symbol)
    name === :light && return Theme()
    name === :dark  || error("presentation_theme: unknown theme :$name (use :light or :dark)")
    return Theme(
        backgroundcolor = :transparent,
        textcolor = :white,
        Axis = (
            backgroundcolor = :transparent,
            xtickcolor       = :white,  ytickcolor       = :white,
            xticklabelcolor  = :white,  yticklabelcolor  = :white,
            xlabelcolor      = :white,  ylabelcolor      = :white,
            titlecolor       = :white,
            bottomspinecolor = :white,  topspinecolor    = :white,
            leftspinecolor   = :white,  rightspinecolor  = :white,
            xgridcolor = (:white, 0.15), ygridcolor = (:white, 0.15),
        ),
        Colorbar = (
            tickcolor   = :white, ticklabelcolor = :white,
            labelcolor  = :white, spinecolor     = :white,
        ),
        Legend = (
            backgroundcolor = :transparent,
            labelcolor = :white, titlecolor = :white,
            framecolor = (:white, 0.4),
        ),
    )
end

# Per-case line colors. Plain primary/secondary names so cases stay
# easy to distinguish on white backgrounds and in printouts. Cycles
# if there are more cases than entries in the palette.
const BASE_CASE_COLORS = [:red, :blue, :green, :orange, :purple,
                           :brown, :magenta, :cyan, :black, :gold]
case_colors = [BASE_CASE_COLORS[mod1(i, length(BASE_CASE_COLORS))]
               for i in 1:length(cases)]

# Standard line styling. Observations are always rendered in light
# grey + dashed + thicker, so the reference is unambiguous against
# any case color but doesn't compete with the model lines for ink.
const CASE_LINEWIDTH = 1.5
const OBS_COLOR      = RGBf(0.55, 0.55, 0.55)
const OBS_LINEWIDTH  = 2.5
const OBS_LINESTYLE  = :dash

run_dir_for(prefix) = "$(prefix)_run"

# ══════════════════════════════════════════════════════════════
# File / time helpers
# ══════════════════════════════════════════════════════════════

function find_first_file(run_dir, prefix, group)
    tag = "$(prefix)_$(group)"
    candidates = filter(f -> startswith(f, tag) && endswith(f, ".jld2") &&
                             !contains(f, "checkpoint"), readdir(run_dir))
    isempty(candidates) && error("No $group files for prefix '$prefix' in $run_dir")
    filename = first(sort(candidates))
    basename_no_part = replace(filename, r"_part\d+" => "")
    return joinpath(run_dir, basename_no_part)
end

function in_window(fts; start_time = 0, stop_time = Inf)
    return findall(t -> start_time <= t <= stop_time, fts.times)
end

function compute_time_mean(fts; start_time = 0, stop_time = Inf)
    idx = in_window(fts; start_time, stop_time)
    isempty(idx) && error("No snapshots in [$start_time, $stop_time]")
    sz  = size(interior(fts[first(idx)]))
    avg = zeros(sz)
    for n in idx
        avg .+= interior(fts[n])
    end
    return avg ./ length(idx)
end

function compute_monthly_means(fts; start_time = 0, stop_time = Inf,
                                reference_date = DateTime(1958, 1, 1))
    idx = in_window(fts; start_time, stop_time)
    isempty(idx) && error("No snapshots in [$start_time, $stop_time]")
    sz     = size(interior(fts[first(idx)]))
    sums   = [zeros(sz) for _ in 1:12]
    counts = zeros(Int, 12)
    for n in idx
        m = month(reference_date + Second(round(Int, fts.times[n])))
        sums[m]  .+= interior(fts[n])
        counts[m] += 1
    end
    return [counts[m] > 0 ? sums[m] ./ counts[m] : nothing for m in 1:12]
end

function compute_mean_and_monthly(fts; start_time = 0, stop_time = Inf,
                                   reference_date = DateTime(1958, 1, 1))
    idx = in_window(fts; start_time, stop_time)
    isempty(idx) && error("No snapshots in [$start_time, $stop_time]")
    sz      = size(interior(fts[first(idx)]))
    total   = zeros(sz)
    monthly = [zeros(sz) for _ in 1:12]
    counts  = zeros(Int, 12)
    for n in idx
        total .+= interior(fts[n])
        m = month(reference_date + Second(round(Int, fts.times[n])))
        monthly[m] .+= interior(fts[n])
        counts[m]   += 1
    end
    mean_out    = total ./ length(idx)
    monthly_out = [counts[m] > 0 ? monthly[m] ./ counts[m] : nothing for m in 1:12]
    return mean_out, monthly_out
end

function cached_download(url; cache_dir = obs_cache_dir,
                              retries = 3,
                              timeout = Inf)
    mkpath(cache_dir)
    path = joinpath(cache_dir, basename(url))
    isfile(path) && return path

    downloader = Downloads.Downloader()
    downloader.easy_hook = (easy, info) ->
        Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_LOW_SPEED_TIME, 0)

    tmp = path * ".part"
    isfile(tmp) && rm(tmp; force = true)

    last_err = nothing
    for attempt in 1:retries
        try
            Downloads.download(url, tmp; timeout, downloader)
            mv(tmp, path; force = true)
            return path
        catch e
            last_err = e
            isfile(tmp) && rm(tmp; force = true)
            if attempt < retries
                delay = 2.0 ^ (attempt - 1)
                @warn "Download failed (attempt $attempt/$retries) — retrying in $(delay)s" url=url error=sprint(showerror, e)
                sleep(delay)
            end
        end
    end
    throw(last_err)
end

# ══════════════════════════════════════════════════════════════
# Grids and masks
# ══════════════════════════════════════════════════════════════

function build_land_mask(grid)
    if grid isa ImmersedBoundaryGrid
        bh = Array(interior(grid.immersed_boundary.bottom_height, :, :, 1))
        return bh .>= 0
    else
        return falses(size(grid, 1), size(grid, 2))
    end
end

function build_ocean_mask_3d(grid)
    Nx, Ny, Nz = size(grid)
    mask = ones(Nx, Ny, Nz)
    if grid isa ImmersedBoundaryGrid
        bh = Array(interior(grid.immersed_boundary.bottom_height, :, :, 1))
        zc = znodes(grid, Center())
        for k in 1:Nz, j in 1:Ny, i in 1:Nx
            zc[k] < bh[i, j] && (mask[i, j, k] = 0.0)
        end
    end
    return mask
end

mask_land!(f, land) = (f[land] .= NaN; f)

# ══════════════════════════════════════════════════════════════
# Plotting helpers
# ══════════════════════════════════════════════════════════════

function panel!(fig, pos, data;
                title="", colormap=:thermal,
                colorrange=nothing, label="",
                nan_color=:lightgray)
    ax = Axis(fig[pos...]; title)
    kw = isnothing(colorrange) ? (;) : (; colorrange)
    hm = heatmap!(ax, data; colormap, nan_color, kw...)
    Colorbar(fig[pos[1], pos[2]+1], hm; label)
    return ax
end

function cpanel!(fig, pos, data;
                 x = nothing, y = nothing,
                 title="", colormap=:thermal,
                 colorrange=nothing, levels=nothing, label="",
                 nan_color=:lightgray, nlevels=21,
                 xlabel="", ylabel="", aspect=nothing)
    axis_kw = isnothing(aspect) ? (; title, xlabel, ylabel) : (; title, xlabel, ylabel, aspect)
    ax = Axis(fig[pos...]; axis_kw...)
    lv = if !isnothing(levels)
        levels
    elseif !isnothing(colorrange)
        range(colorrange[1], colorrange[2]; length=nlevels)
    else
        nothing
    end
    kw = isnothing(lv) ? (;) : (; levels=lv)
    xs = isnothing(x) ? (1:size(data, 1)) : x
    ys = isnothing(y) ? (1:size(data, 2)) : y
    hm = contourf!(ax, xs, ys, data;
                   colormap, nan_color,
                   extendlow=:auto, extendhigh=:auto, kw...)
    Colorbar(fig[pos[1], pos[2]+1], hm; label)
    return ax
end

# ── GeoMakie contour-filled panel ─────────────────────────────────────
#
# `geo_panel!` is the geographic-projection counterpart to `cpanel!`.
# It builds a `GeoAxis` with the requested PROJ string, overlays
# Natural-Earth coastlines + land, and renders the data as filled
# contours with `extendlow/high = :auto` so values outside the level
# range still get the end-of-colormap color (no white margins at the
# extremes).
#
# `x` and `y` are physical longitudes / latitudes (degrees). `data` is
# `(Nx, Ny)`. NaN cells render transparent so the underlying land /
# coastline shows through naturally — no `nan_color` machinery needed.

# Shift a 0..360-longitude axis + data to -180..180 (sorted ascending).
# GeoMakie's projections (e.g. Robinson) clip / silently drop data points
# with longitudes outside [-180, 180]. Our regridded fields use the
# Oceananigans 0..360 convention, so we wrap before passing to contourf.
function to_minus180_180(lon, data)
    any(lon .> 180) || return lon, data
    split = findfirst(>(180), lon)
    new_lon  = vcat(lon[split:end] .- 360, lon[1:split-1])
    new_data = vcat(data[split:end, :],   data[1:split-1, :])
    return new_lon, new_data
end

# Sutherland-Hodgman polygon clipping against a horizontal lat line.
# `keep_above = true`  → drop points with lat < lat_threshold;
# `keep_above = false` → drop points with lat > lat_threshold.
function sh_clip_one_side(pts, lat_threshold, keep_above::Bool)
    isempty(pts) && return pts
    is_inside(p) = keep_above ? (p[2] >= lat_threshold) : (p[2] <= lat_threshold)
    out = similar(pts, 0)
    n = length(pts)
    for i in 1:n
        cur = pts[i]
        prev = pts[mod1(i-1, n)]
        cur_in  = is_inside(cur)
        prev_in = is_inside(prev)
        if cur_in
            if !prev_in
                t = (lat_threshold - prev[2]) / (cur[2] - prev[2])
                push!(out, typeof(cur)(prev[1] + t*(cur[1]-prev[1]), lat_threshold))
            end
            push!(out, cur)
        elseif prev_in
            t = (lat_threshold - prev[2]) / (cur[2] - prev[2])
            push!(out, typeof(cur)(prev[1] + t*(cur[1]-prev[1]), lat_threshold))
        end
    end
    return out
end

# Clip a GeometryBasics polygon to [lat_min, lat_max]. Returns
# `nothing` if the result degenerates to fewer than 3 vertices.
function clip_polygon_to_lat_band(polygon, lat_min, lat_max)
    pts = collect(GeometryBasics.coordinates(polygon))
    pts = sh_clip_one_side(pts, lat_max, false)   # drop lat > lat_max
    pts = sh_clip_one_side(pts, lat_min, true)    # drop lat < lat_min
    length(pts) < 3 && return nothing
    return GeometryBasics.Polygon(pts)
end

# Cached lat-band-clipped Natural-Earth land. Keyed by (lat_min, lat_max).
# Polar panels otherwise fill their rectangular-axis corners with the
# (projected) low-latitude portion of the polygon — Europe/N.America in
# NH, sparse in SH — producing a NH-vs-SH grey/white asymmetry. Clipping
# the polygon to the panel's lat band eliminates that contribution while
# preserving the in-disk land fill.
const CLIPPED_LAND_CACHE = Dict{Tuple{Float64, Float64}, Any}()

function clipped_land(lat_min, lat_max)
    key = (Float64(lat_min), Float64(lat_max))
    return get!(CLIPPED_LAND_CACHE, key) do
        polys = GeoMakie.land()
        out = eltype(polys)[]
        for p in polys
            cp = clip_polygon_to_lat_band(p, lat_min, lat_max)
            cp === nothing || push!(out, cp)
        end
        out
    end
end

function geo_panel!(fig, pos, data;
                    x, y,
                    projection = "+proj=robin",
                    title = "", colormap = :thermal,
                    colorrange = nothing, levels = nothing, label = "",
                    nlevels = 21,
                    coastlines = true,
                    coast_color = :black, coast_linewidth = 0.4,
                    land_color = :lightgray,
                    lonlims = nothing, latlims = nothing,
                    obs_contour = nothing,
                    obs_levels = [0.15],
                    obs_color = :black,
                    obs_linewidth = 1.5)
    x_in = x
    x, data = to_minus180_180(x, data)
    # `titlegap` clears the graticule labels at the disk edge — without
    # it the title text in polar-stereographic panels can overlap the
    # longitude labels (e.g. "30°") at the top.
    ga = GeoAxis(fig[pos...]; dest = projection, title, titlegap = 24)
    # GeoMakie's xlims!/ylims! take longitude/latitude in degrees, NOT
    # projected coordinates. Used to clip polar caps so the data fills
    # the panel instead of being squashed into a whole-globe view.
    isnothing(lonlims) || xlims!(ga, lonlims...)
    isnothing(latlims) || ylims!(ga, latlims...)
    # Land polygon is drawn UNDER the contourf so NaN cells (regridded
    # land) show the polygon. Data is drawn on top — it should be
    # NaN-masked at land cells; if it isn't, the bleed is a regridder
    # bug, not a rendering one. (Drawing the polygon on top doesn't
    # work in polar-stereographic projections — the Natural-Earth
    # polygon edges wrap incorrectly and fill the whole disk.)
    #
    # When `latlims` is set (polar panels), pre-clip the polygon to
    # the latitude band so the low-latitude land (Europe / N. America)
    # doesn't fill the rectangular axis corners outside the disk.
    if !isnothing(land_color)
        land_polys = isnothing(latlims) ? GeoMakie.land() : clipped_land(latlims...)
        poly!(ga, land_polys; color = land_color, strokewidth = 0)
    end
    lv = if !isnothing(levels)
        levels
    elseif !isnothing(colorrange)
        range(colorrange[1], colorrange[2]; length = nlevels)
    else
        nothing
    end
    kw = isnothing(lv) ? (;) : (; levels = lv)
    hm = contourf!(ga, x, y, data;
                   colormap,
                   extendlow = :auto, extendhigh = :auto, kw...)
    if coastlines
        lines!(ga, GeoMakie.coastlines(); color = coast_color,
               linewidth = coast_linewidth)
    end
    # Obs reference contour (e.g. 15% SIC edge) drawn last so it sits
    # on top of the filled data and the coastlines.
    if !isnothing(obs_contour)
        _, obs_shifted = to_minus180_180(x_in, obs_contour)
        contour!(ga, x, y, obs_shifted;
                 levels = obs_levels, color = obs_color, linewidth = obs_linewidth)
    end
    Colorbar(fig[pos[1], pos[2] + 1], hm; label)
    return ga
end

# ══════════════════════════════════════════════════════════════
# ECCO4 free-surface climatology
# ══════════════════════════════════════════════════════════════

function ecco_ssh_climatology_native(; start_date = DateTime(1992, 1, 1),
                                       end_date   = DateTime(2012, 12, 1),
                                       cache_dir  = obs_cache_dir)
    cache_file = joinpath(cache_dir, "ssh_ecco4_$(year(start_date))_$(year(end_date))_native.jld2")
    dates = start_date:Month(1):end_date
    if isfile(cache_file)
        return JLD2.load(cache_file, "ssh_mean")
    end
    @info "  Computing ECCO4 SSH climatology over $(length(dates)) months (one-time)..."
    first_field = Field(Metadatum(:free_surface; dataset = ECCO4Monthly(), date = first(dates)), CPU())
    ssh_mean    = copy(Array(interior(first_field)))
    for date in dates[2:end]
        f = Field(Metadatum(:free_surface; dataset = ECCO4Monthly(), date), CPU())
        ssh_mean .+= Array(interior(f))
    end
    ssh_mean ./= length(dates)
    JLD2.jldsave(cache_file; ssh_mean)
    return ssh_mean
end

const ECCO_SSH_NATIVE_MEAN_REF = Ref{Any}(nothing)

function ecco_ssh_on_grid(grid; reference_date = DateTime(1992, 1, 1))
    if isnothing(ECCO_SSH_NATIVE_MEAN_REF[])
        ECCO_SSH_NATIVE_MEAN_REF[] = ecco_ssh_climatology_native()
    end
    template = Field(Metadatum(:free_surface; dataset = ECCO4Monthly(), date = reference_date), CPU())
    interior(template) .= ECCO_SSH_NATIVE_MEAN_REF[]
    dst = Field{Center, Center, Nothing}(grid)
    interpolate!(dst, template)
    return dropdims(Array(interior(dst)); dims = 3)
end

# ══════════════════════════════════════════════════════════════
# de Boyer Montégut MLD climatology
# ══════════════════════════════════════════════════════════════

const DBM_MLD_URL = get(ENV, "DBM_MLD_URL", "https://mld.ifremer.fr/data/mld_DR003_c1m_reg2.0.nc")

function centers_to_edges(centers; clamp_to = nothing)
    Δfirst = centers[2] - centers[1]
    Δlast  = centers[end] - centers[end-1]
    edges  = Vector{Float64}(undef, length(centers) + 1)
    edges[1]   = centers[1] - Δfirst / 2
    edges[end] = centers[end] + Δlast / 2
    for i in 2:length(centers)
        edges[i] = (centers[i-1] + centers[i]) / 2
    end
    if !isnothing(clamp_to)
        lo, hi = clamp_to
        edges[1]   = max(edges[1], lo)
        edges[end] = min(edges[end], hi)
    end
    return edges
end

function dbm_mld_climatology_on_grid(grid;
                                     file = get(ENV, "DBM_MLD_FILE", joinpath(obs_cache_dir, basename(DBM_MLD_URL))),
                                     var  = get(ENV, "DBM_MLD_VAR", "mld"))
    if !isfile(file)
        try
            @info "  Downloading dBM MLD climatology from $DBM_MLD_URL"
            file = cached_download(DBM_MLD_URL)
        catch e
            @warn "dBM MLD auto-download failed — skipping reference. Manually download from https://mld.ifremer.fr/Surface_Mixed_Layer_Depth.php and set DBM_MLD_FILE." error=sprint(showerror, e)
            return nothing
        end
    end
    ds = NCDatasets.NCDataset(file)
    mld_raw = Array(ds[var][:, :, :])
    lon_vec = Float64.(Array(ds["lon"][:]))
    lat_vec = Float64.(Array(ds["lat"][:]))
    close(ds)

    mld_raw = Float64.(coalesce.(mld_raw, NaN))
    mld_raw[mld_raw .> 1e8] .= NaN

    lon_edges = centers_to_edges(lon_vec)
    lat_edges = centers_to_edges(lat_vec; clamp_to = (-90, 90))
    Nlon, Nlat, Nm = size(mld_raw)

    src_grid = LatitudeLongitudeGrid(CPU();
        size = (Nlon, Nlat, 1),
        longitude = lon_edges,
        latitude  = lat_edges,
        z = (0, 1))

    src = Field{Center, Center, Nothing}(src_grid)
    dst = Field{Center, Center, Nothing}(grid)
    Nx, Ny = size(grid, 1), size(grid, 2)
    out = Array{Float64, 3}(undef, Nx, Ny, Nm)
    for m in 1:Nm
        clean = replace(mld_raw[:, :, m], NaN => 0.0)
        interior(src) .= reshape(clean, Nlon, Nlat, 1)
        interpolate!(dst, src)
        out[:, :, m] = Array(interior(dst))[:, :, 1]
    end
    return out
end

# ══════════════════════════════════════════════════════════════
# NCEP/NCAR Reanalysis 1 wind-stress climatology
# ══════════════════════════════════════════════════════════════

const NCEP_TAUU_URL = get(ENV, "NCEP_TAUU_URL", "https://psl.noaa.gov/thredds/fileServer/Datasets/ncep.reanalysis.derived/surface_gauss/uflx.sfc.mon.ltm.nc")
const NCEP_TAUV_URL = get(ENV, "NCEP_TAUV_URL", "https://psl.noaa.gov/thredds/fileServer/Datasets/ncep.reanalysis.derived/surface_gauss/vflx.sfc.mon.ltm.nc")

function ncep_wind_stress_on_grid(grid;
                                  tauu_file = get(ENV, "NCEP_TAUU_FILE", joinpath(obs_cache_dir, basename(NCEP_TAUU_URL))),
                                  tauv_file = get(ENV, "NCEP_TAUV_FILE", joinpath(obs_cache_dir, basename(NCEP_TAUV_URL))),
                                  tauu_var  = get(ENV, "NCEP_TAUU_VAR",  "uflx"),
                                  tauv_var  = get(ENV, "NCEP_TAUV_VAR",  "vflx"))
    function ensure(file, url)
        isfile(file) && return file
        @info "  Downloading NCEP wind-stress climatology from $url"
        return cached_download(url)
    end
    try
        tauu_file = ensure(tauu_file, NCEP_TAUU_URL)
        tauv_file = ensure(tauv_file, NCEP_TAUV_URL)
    catch e
        @warn "NCEP auto-download failed — skipping wind-stress reference. Provide netCDFs via NCEP_TAUU_FILE / NCEP_TAUV_FILE." error=sprint(showerror, e)
        return nothing, nothing
    end

    function read_stress(file, var)
        ds = NCDatasets.NCDataset(file)
        raw = Array(ds[var])
        lon_vec = Float64.(Array(ds["lon"][:]))
        lat_vec = Float64.(Array(ds["lat"][:]))
        close(ds)
        for d in reverse(findall(==(1), size(raw)))
            raw = dropdims(raw; dims=d)
        end
        return Float64.(coalesce.(raw, NaN)), lon_vec, lat_vec
    end

    τu_raw, lon_vec, lat_vec = read_stress(tauu_file, tauu_var)
    τv_raw, _, _             = read_stress(tauv_file, tauv_var)

    annual_mean(f) = ndims(f) == 3 ? dropdims(mean(f, dims=3), dims=3) : f
    τx_2d = -annual_mean(τu_raw)
    τy_2d = -annual_mean(τv_raw)

    if lat_vec[1] > lat_vec[end]
        lat_vec = reverse(lat_vec)
        τx_2d   = reverse(τx_2d; dims=2)
        τy_2d   = reverse(τy_2d; dims=2)
    end

    lon_edges = centers_to_edges(lon_vec)
    lat_edges = centers_to_edges(lat_vec; clamp_to = (-90, 90))
    Nlon, Nlat = length(lon_vec), length(lat_vec)

    src_grid = LatitudeLongitudeGrid(CPU();
        size = (Nlon, Nlat, 1),
        longitude = lon_edges,
        latitude  = lat_edges,
        z = (0, 1))

    src = Field{Center, Center, Nothing}(src_grid)
    dst = Field{Center, Center, Nothing}(grid)
    function regrid(data_2d)
        clean = replace(data_2d, NaN => 0.0)
        interior(src) .= reshape(clean, Nlon, Nlat, 1)
        interpolate!(dst, src)
        return dropdims(Array(interior(dst)); dims=3)
    end
    return regrid(τx_2d), regrid(τy_2d)
end

# ══════════════════════════════════════════════════════════════
# Postprocess MLD with a configurable reference depth (vertical interp only)
# ══════════════════════════════════════════════════════════════

"""
    mld_with_reference_depth(buoyancy, grid;
                             ocean_mask = nothing,
                             reference_depth = 10.0,
                             Δb★ = 2.87e-4)

Compute mixed-layer depth at every (i, j) by vertically interpolating
the buoyancy field at `z = -reference_depth` and walking down the
column until the buoyancy drop relative to that reference exceeds
`Δb★` (de Boyer Montégut DR003 default at ρ₀ = 1025 kg/m³). Returns a
2-D array of MLD in meters; NaN where the column is land.

`ocean_mask` is a 3-D `(Nx, Ny, Nz)` array with `1` for ocean and `0`
for land cells (e.g. as returned by `build_ocean_mask_3d(grid)`). When
omitted the routine falls back to NaN-only land detection, which can
mistake a legitimately zero buoyancy in the open ocean for land — pass
the mask whenever it is available.

Pure column-wise computation — no horizontal interpolation needed
because the z-coordinate is the same at every (i, j).
"""
function mld_with_reference_depth(buoyancy::AbstractArray{<:Any, 3}, grid;
                                  ocean_mask = nothing,
                                  reference_depth = 10.0,
                                  Δb★ = 2.87e-4)
    Nx, Ny, Nz = size(buoyancy)
    z   = collect(znodes(grid, Center()))
    zʳ  = -reference_depth
    mld = fill(NaN, Nx, Ny)

    is_land(i, j, k) = isnothing(ocean_mask) ? false : ocean_mask[i, j, k] == 0

    @inbounds for j in 1:Ny, i in 1:Nx
        is_land(i, j, Nz) && continue
        bN = buoyancy[i, j, Nz]
        isnan(bN) && continue

        # Bracket zʳ between cell centers k⁻ (deeper) and k⁺ (shallower).
        k⁺ = Nz
        while k⁺ > 1 && z[k⁺ - 1] >= zʳ
            k⁺ -= 1
        end
        bʳ = if k⁺ == 1
            bN
        else
            k⁻ = k⁺ - 1
            w  = (zʳ - z[k⁻]) / (z[k⁺] - z[k⁻])
            buoyancy[i, j, k⁻] * (1 - w) + buoyancy[i, j, k⁺] * w
        end

        # Walk down from zʳ; linearly interpolate where Δb crosses Δb★.
        zₚ  = zʳ
        Δbₚ = zero(eltype(buoyancy))
        z★  = NaN
        for k in (Nz - 1):-1:1
            z[k] >= zʳ && continue
            is_land(i, j, k) && break
            bk = buoyancy[i, j, k]
            isnan(bk) && break
            Δb = bʳ - bk
            if Δb >= Δb★
                z★ = zₚ + (Δb★ - Δbₚ) / (Δb - Δbₚ) * (z[k] - zₚ)
                break
            end
            zₚ = z[k]; Δbₚ = Δb
        end

        mld[i, j] = isnan(z★) ? -zₚ : -z★
    end

    return mld
end

# ══════════════════════════════════════════════════════════════
# Sea-ice integral helpers (used by the :ice_diag loader)
# ══════════════════════════════════════════════════════════════

arctic_condition(i, j, k, grid, args...)    = φnode(i, j, k, grid, Center(), Center(), Center()) > 0
antarctic_condition(i, j, k, grid, args...) = φnode(i, j, k, grid, Center(), Center(), Center()) < 0

function compute_ice_diagnostics(run_dir, prefix, grid;
                                 start_time = 0, stop_time = Inf,
                                 reference_date = DateTime(1958, 1, 1),
                                 extent_threshold = 0.15)
    surface_file      = find_first_file(run_dir, prefix, "surface")
    thickness_fts     = FieldTimeSeries(surface_file, "sithick"; backend = deepcopy(FTS_BACKEND))
    concentration_fts = FieldTimeSeries(surface_file, "siconc";  backend = deepcopy(FTS_BACKEND))

    # Restrict the per-snapshot integration loop to the case averaging
    # window — pre-window snapshots can't enter any monthly bin so
    # there's no reason to integrate them.
    idx = findall(t -> start_time <= t <= stop_time, thickness_fts.times)
    isempty(idx) && error("compute_ice_diagnostics: no snapshots in [$start_time, $stop_time]")
    Nidx = length(idx)
    @info "  $prefix: integrating sea-ice over $Nidx snapshots..."

    snapshot_dates   = [reference_date + Second(round(Int, thickness_fts.times[n])) for n in idx]
    arctic_volume    = zeros(Nidx)
    antarctic_volume = zeros(Nidx)
    arctic_extent    = zeros(Nidx)
    antarctic_extent = zeros(Nidx)
    arctic_area      = zeros(Nidx)
    antarctic_area   = zeros(Nidx)

    thickness     = Field{Center, Center, Nothing}(grid)
    concentration = Field{Center, Center, Nothing}(grid)
    ice_volume    = Field{Center, Center, Nothing}(grid)
    extent_mask   = Field{Center, Center, Nothing}(grid)

    arctic_vol_int     = Field(Integral(ice_volume;    condition = arctic_condition))
    antarctic_vol_int  = Field(Integral(ice_volume;    condition = antarctic_condition))
    arctic_area_int    = Field(Integral(concentration; condition = arctic_condition))
    antarctic_area_int = Field(Integral(concentration; condition = antarctic_condition))
    arctic_ext_int     = Field(Integral(extent_mask;   condition = arctic_condition))
    antarctic_ext_int  = Field(Integral(extent_mask;   condition = antarctic_condition))

    for (k, n) in enumerate(idx)
        set!(thickness,     thickness_fts[n])
        set!(concentration, concentration_fts[n])
        interior(ice_volume) .= interior(thickness) .* interior(concentration)

        compute!(arctic_vol_int);  compute!(antarctic_vol_int)
        arctic_volume[k]    = arctic_vol_int[1, 1, 1]
        antarctic_volume[k] = antarctic_vol_int[1, 1, 1]

        compute!(arctic_area_int); compute!(antarctic_area_int)
        arctic_area[k]    = arctic_area_int[1, 1, 1]
        antarctic_area[k] = antarctic_area_int[1, 1, 1]

        # Build the > threshold extent mask in place — no per-snapshot
        # temporary array.
        interior(extent_mask) .= interior(concentration) .> extent_threshold
        compute!(arctic_ext_int); compute!(antarctic_ext_int)
        arctic_extent[k]    = arctic_ext_int[1, 1, 1]
        antarctic_extent[k] = antarctic_ext_int[1, 1, 1]
    end

    months_used = month.(snapshot_dates)
    monthly(field) = [mean(field[months_used .== m]) for m in 1:12]

    return (; arctic_volume, antarctic_volume,
              arctic_extent, antarctic_extent,
              arctic_area, antarctic_area, snapshot_dates,
              arctic_volume_monthly    = monthly(arctic_volume),
              antarctic_volume_monthly = monthly(antarctic_volume),
              arctic_extent_monthly    = monthly(arctic_extent),
              antarctic_extent_monthly = monthly(antarctic_extent),
              arctic_area_monthly      = monthly(arctic_area),
              antarctic_area_monthly   = monthly(antarctic_area))
end

# ══════════════════════════════════════════════════════════════
# Observational sea-ice climatologies (NSIDC, PIOMAS) — global cache.
# ══════════════════════════════════════════════════════════════

function load_piomas_monthly()
    url   = "https://psc.apl.uw.edu/wordpress/wp-content/uploads/schweiger/ice_volume/PIOMAS.monthly.Current.v2.1.csv"
    raw   = readdlm(cached_download(url), ','; skipstart=1)
    vol   = Float64.(raw[:, 2:13])
    vol[vol .== -1] .= NaN
    volume_monthly     = vec(mapslices(x -> mean(filter(!isnan, x)), vol; dims=1))
    volume_monthly_std = vec(mapslices(vol; dims=1) do x
        y = filter(!isnan, x)
        length(y) > 1 ? std(y) : 0.0
    end)
    return (; volume_monthly, volume_monthly_std)
end

function load_nsidc(hemisphere)
    prefix = hemisphere == "north" ? "N" : "S"
    extent_monthly     = zeros(12)
    extent_monthly_std = zeros(12)
    area_monthly       = zeros(12)
    area_monthly_std   = zeros(12)
    for m in 1:12
        url = "https://noaadata.apps.nsidc.org/NOAA/G02135/$(hemisphere)/monthly/data/$(prefix)_$(lpad(m, 2, '0'))_extent_v4.0.csv"
        raw = readlines(cached_download(url))
        extents = Float64[]; areas = Float64[]
        for line in raw
            parts = split(line, ',')
            length(parts) >= 6 || continue
            ext = tryparse(Float64, strip(parts[5]))
            ar  = tryparse(Float64, strip(parts[6]))
            (isnothing(ext) || ext == -9999) && continue
            (isnothing(ar)  || ar  == -9999) && continue
            push!(extents, ext); push!(areas, ar)
        end
        extent_monthly[m]     = mean(extents)
        extent_monthly_std[m] = length(extents) > 1 ? std(extents) : 0.0
        area_monthly[m]       = mean(areas)
        area_monthly_std[m]   = length(areas) > 1 ? std(areas) : 0.0
    end
    return (; extent_monthly, extent_monthly_std, area_monthly, area_monthly_std)
end

# Global per-process cache for observational climatologies. Each
# accessor returns `nothing` if the download fails (mirrors the dBM /
# NCEP convention) so figures degrade gracefully to a model-only plot.
# Sentinel `:download_failed` distinguishes "download failed once,
# don't retry this session" from "not yet attempted".
const NSIDC_NORTH_REF = Ref{Any}(nothing)
const NSIDC_SOUTH_REF = Ref{Any}(nothing)
const PIOMAS_REF      = Ref{Any}(nothing)

function try_load(ref, label, builder)
    ref[] === :download_failed && return nothing
    isnothing(ref[]) || return ref[]
    try
        ref[] = builder()
    catch err
        @warn "$label download failed — skipping reference line." error = sprint(showerror, err)
        ref[] = :download_failed
        return nothing
    end
    return ref[]
end

nsidc_arctic()    = try_load(NSIDC_NORTH_REF, "NSIDC (north)", () -> load_nsidc("north"))
nsidc_antarctic() = try_load(NSIDC_SOUTH_REF, "NSIDC (south)", () -> load_nsidc("south"))
piomas_monthly()  = try_load(PIOMAS_REF,      "PIOMAS",        load_piomas_monthly)

# ══════════════════════════════════════════════════════════════
# HadISST1 sea-ice concentration climatology (Met Office Hadley Centre)
# ══════════════════════════════════════════════════════════════
#
# HadISST1 (Rayner et al. 2003) blends satellite passive-microwave SIC
# with earlier sources back to 1870 on a global 1°×1° lat-lon grid.
# Distributed as a single gzipped NetCDF (`HadISST_ice.nc.gz`, ~30 MB);
# we gunzip on first download into the obs cache and compute the
# monthly climatology over a fixed window (1979–2007 to match the
# Adcroft 2019 OM4 reference window).
#
# Native HadISST conventions: longitude −179.5..179.5 (ascending),
# latitude 89.5..−89.5 (descending), `sic` as 0..1 fraction with land /
# missing flagged at −1 / −1e30. Because the resolution matches the
# shared `LATLON_*_CENTERS` 1° grid exactly, we line them up by axis
# remap (longitude roll + latitude flip) — no interpolation.

const HADISST_SIC_URL = get(ENV, "HADISST_SIC_URL",
    "https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_ice.nc.gz")

function gunzip_to_sibling(path_gz::AbstractString)
    endswith(path_gz, ".gz") || return path_gz
    out = path_gz[1:end-3]
    isfile(out) && return out
    @info "  Gunzipping $(basename(path_gz))"
    # `gunzip -k` keeps the .gz file but isn't supported by older gzip
    # builds (e.g. some cluster login nodes). Use `gunzip -c` (write to
    # stdout) with a Julia pipeline redirect — works on every gzip.
    run(pipeline(`gunzip -c $path_gz`, stdout = out))
    return out
end

function load_hadisst_sic_climatology(; start_year = 1979, end_year = 2007,
                                       cache_dir = obs_cache_dir)
    cache_file = joinpath(cache_dir, "sic_hadisst_$(start_year)_$(end_year).jld2")
    if isfile(cache_file)
        return JLD2.load(cache_file, "monthly")
    end

    file = gunzip_to_sibling(cached_download(HADISST_SIC_URL))
    @info "  Reading HadISST monthly SIC and building $(start_year)–$(end_year) climatology..."

    ds       = NCDatasets.NCDataset(file)
    lon_vec  = Float64.(Array(ds["longitude"][:]))
    lat_vec  = Float64.(Array(ds["latitude"][:]))
    time_vec = Array(ds["time"][:])
    Nlon     = length(lon_vec)
    Nlat     = length(lat_vec)
    sic_var  = ds["sic"]
    yrs      = Dates.year.(time_vec)
    mns      = Dates.month.(time_vec)
    keep     = findall(y -> start_year <= y <= end_year, yrs)
    isempty(keep) && (close(ds); error("HadISST: no months in [$start_year, $end_year]"))

    sums      = [zeros(Nlon, Nlat) for _ in 1:12]
    good_cnt  = [zeros(Int, Nlon, Nlat) for _ in 1:12]
    auto_scale_pending = true   # decide 0..1 vs 0..100 on first finite slab

    for n in keep
        slab = Float64.(coalesce.(Array(sic_var[:, :, n]), NaN))
        # Land / fill values
        @inbounds for j in 1:Nlat, i in 1:Nlon
            v = slab[i, j]
            (v < -0.5 || v > 200.0) && (slab[i, j] = NaN)
        end
        if auto_scale_pending
            finite_max = maximum(Iterators.filter(isfinite, slab); init = 0.0)
            if finite_max > 1.5
                slab ./= 100
            end
            auto_scale_pending = false
        end
        m = mns[n]
        @inbounds for j in 1:Nlat, i in 1:Nlon
            v = slab[i, j]
            if isfinite(v)
                sums[m][i, j]     += v
                good_cnt[m][i, j] += 1
            end
        end
    end
    close(ds)

    # Native-grid monthly mean (NaN where every contributing month was NaN).
    monthly_native = Vector{Matrix{Float64}}(undef, 12)
    for m in 1:12
        out = fill(NaN, Nlon, Nlat)
        @inbounds for j in 1:Nlat, i in 1:Nlon
            good_cnt[m][i, j] > 0 && (out[i, j] = sums[m][i, j] / good_cnt[m][i, j])
        end
        monthly_native[m] = out
    end

    # Remap HadISST axes (-180..180, lat descending) to the shared 1° lat-lon
    # convention (0..360, lat ascending). HadISST native grid is exactly 1°×1°
    # — same as `LATLON_*_CENTERS` — so this is a pure sort, no interpolation.
    @assert Nlon == length(LATLON_LON_CENTERS) "HadISST longitude grid mismatch: $Nlon vs $(length(LATLON_LON_CENTERS))"
    @assert Nlat == length(LATLON_LAT_CENTERS) "HadISST latitude grid mismatch: $Nlat vs $(length(LATLON_LAT_CENTERS))"
    i_perm = sortperm(mod.(lon_vec, 360))
    j_perm = sortperm(lat_vec)

    monthly = [monthly_native[m][i_perm, j_perm] for m in 1:12]
    JLD2.jldsave(cache_file; monthly)
    return monthly
end

const HADISST_SIC_REF = Ref{Any}(nothing)

hadisst_sic_climatology() =
    try_load(HADISST_SIC_REF, "HadISST SIC climatology", load_hadisst_sic_climatology)

# ══════════════════════════════════════════════════════════════
# RAPID-MOCHA 26.5°N AMOC climatology + time series
# ══════════════════════════════════════════════════════════════
#
# The RAPID-MOCHA array (McCarthy et al. 2015; Moat et al. 2023) has
# monitored the AMOC at 26.5°N continuously since April 2004. We use
# two NetCDFs from the project:
#
#   * `moc_vertical.nc`  — daily vertical streamfunction ψ(z, t) (Sv)
#                          on a fixed depth axis (positive m).
#   * `moc_transports.nc` — daily AMOC components; `moc_mar_hc10` is
#                           the upper-cell maximum, used as the AMOC
#                           index time series.
#
# The versionless `rapid_data/` path is the canonical home for the
# current release (v2024.1a as of 2026-01); rapid.ac.uk has previously
# served the same files from snapshot-dated `YYYY-MM/` folders that
# disappear when a new release is cut. Override with `RAPID_*_URL`
# env vars if the path ever moves again.

const RAPID_VERTICAL_URL = get(ENV, "RAPID_VERTICAL_URL",
    "https://rapid.ac.uk/sites/default/files/rapid_data/moc_vertical.nc")
const RAPID_TRANSPORTS_URL = get(ENV, "RAPID_TRANSPORTS_URL",
    "https://rapid.ac.uk/sites/default/files/rapid_data/moc_transports.nc")

function nanmean(v)
    finite = filter(isfinite, v)
    return isempty(finite) ? NaN : mean(finite)
end
function nanstd(v)
    finite = filter(isfinite, v)
    return length(finite) < 2 ? NaN : std(finite)
end

# Climatological ψ(z) at 26.5°N: time-mean and standard deviation of
# the daily RAPID streamfunction profile. NaN at depths/times where
# the array did not measure.
function load_rapid_amoc_profile(; cache_dir = obs_cache_dir)
    file = cached_download(RAPID_VERTICAL_URL; cache_dir)
    @info "  Reading RAPID vertical streamfunction at 26.5°N..."
    ds = NCDatasets.NCDataset(file)
    depth_var = haskey(ds, "depth") ? ds["depth"] : ds["z"]
    psi_var   = haskey(ds, "stream_function_mar") ? ds["stream_function_mar"] :
                haskey(ds, "stream_function") ? ds["stream_function"] :
                error("RAPID moc_vertical.nc: cannot find streamfunction variable")
    depth = Float64.(Array(depth_var[:]))
    raw   = Float64.(coalesce.(Array(psi_var[:, :]), NaN))   # in Sv
    # The v2024.1a release declares `stream_function_mar(time, depth)`,
    # which NCDatasets returns as a Julia array sized (Ntime, Ndepth);
    # older releases sometimes flip the order. Orient on the depth axis
    # so `psi[k, :]` is unambiguously the time series at depth k.
    psi = size(raw, 1) == length(depth) ? raw : permutedims(raw, (2, 1))
    close(ds)
    mean_psi = [nanmean(@view psi[k, :]) for k in 1:length(depth)]
    std_psi  = [nanstd(@view psi[k, :])  for k in 1:length(depth)]
    return (depth = depth, psi_mean = mean_psi, psi_std = std_psi)
end

# Monthly-binned AMOC max time series from `moc_mar_hc10`. Returned
# coordinate is decimal calendar year so it can be plotted alongside a
# model time series that has been converted to year-since-1958.
function load_rapid_amoc_timeseries(; cache_dir = obs_cache_dir)
    file = cached_download(RAPID_TRANSPORTS_URL; cache_dir)
    @info "  Reading RAPID AMOC index time series..."
    ds = NCDatasets.NCDataset(file)
    time_var = Array(ds["time"][:])
    moc_name = haskey(ds, "moc_mar_hc10") ? "moc_mar_hc10" :
               haskey(ds, "moc_mar_hc")   ? "moc_mar_hc"   :
               error("RAPID moc_transports.nc: cannot find AMOC index variable")
    psi_max  = Float64.(coalesce.(Array(ds[moc_name][:]), NaN))
    close(ds)
    yrs = Dates.year.(time_var)
    mns = Dates.month.(time_var)
    keys_ym = sort!(unique(collect(zip(yrs, mns))))
    centers = Float64[]
    means   = Float64[]
    stds    = Float64[]
    for (y, m) in keys_ym
        idx = findall(i -> yrs[i] == y && mns[i] == m, eachindex(time_var))
        vals = filter(isfinite, view(psi_max, idx))
        isempty(vals) && continue
        push!(centers, y + (m - 0.5) / 12)
        push!(means, mean(vals))
        push!(stds, length(vals) > 1 ? std(vals) : 0.0)
    end
    return (year = centers, psi_max = means, psi_max_std = stds)
end

const RAPID_PROFILE_REF    = Ref{Any}(nothing)
const RAPID_TIMESERIES_REF = Ref{Any}(nothing)

rapid_amoc_profile()    = try_load(RAPID_PROFILE_REF,    "RAPID ψ(z) at 26.5°N", load_rapid_amoc_profile)
rapid_amoc_timeseries() = try_load(RAPID_TIMESERIES_REF, "RAPID AMOC index",     load_rapid_amoc_timeseries)
