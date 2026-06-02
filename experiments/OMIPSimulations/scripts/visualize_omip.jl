#!/usr/bin/env julia
# visualize_omip.jl — OMIP diagnostic figures: set up + render.
#
# This file is designed for two workflows.
#
# REPL (preferred for iteration)
# ------------------------------
#     julia> include("visualize_omip.jl")
#     # nothing renders; `caches`, `labels`, `cases` and all fig01..fig21
#     # are now defined.
#     julia> fig04(caches, labels, cases)        # render just fig 4
#     julia> fig17(caches, labels, cases)        # later, render fig 17
#     julia> fig04(caches, labels, cases)        # rerun: cache hits, fast
#
# Script (batch)
# --------------
#     julia visualize_omip.jl                    # render all 21 figs
#     FIG=4 julia visualize_omip.jl              # render only fig 4
#     FIG=1,4,7 julia visualize_omip.jl          # comma list
#     FIG=14-19 julia visualize_omip.jl          # range
#     FIG=1,3,9-12,21 julia visualize_omip.jl    # mix
#     THEME=dark julia visualize_omip.jl         # transparent bg, white axes
#     julia visualize_omip.jl my_output_dir      # custom output dir
#
# How sharing works
# -----------------
# Each figure declares its data needs implicitly by calling
# `get_field(cache, :sym)`. Loaders form a DAG (e.g. `:sst_bias` ←
# `:sst` ← `:tos_fts`, plus `:woa_temperature`). Per orchestrator
# session, every loader fires at most once per case — so running
# `fig01` then `fig02` reads `tos` and `sos` once each but never
# reloads the WOA file.
#
# Running `fig01(caches, labels, cases)` alone touches only
# `:tos_fts` + `:woa_temperature` + the masks/grid; nothing else is
# loaded. That's the
# "minimum-time, isolation" property the refactor was designed for.
#
# Edit the `cases` list below before the first include.
#
# Each case specifies its averaging window in one of two mutually
# exclusive ways:
#   - absolute: `start_time` and/or `stop_time` in seconds (omit either
#     to default to 0 / Inf).
#   - relative: `years_from_end = N` → average the last N years of the
#     run, computed from the latest snapshot in the surface JLD2 file.

# ══════════════════════════════════════════════════════════════
# Configuration
# ══════════════════════════════════════════════════════════════

cases = [
    (prefix = "orca_corrected_snow_rbvd_bih50days",             label = "ORCA RBVD",           years_from_end = 5),
    # (prefix = "orca_ncar_snow",                                 label = "ORCA New NCAR",       years_from_end = 5),
    # (prefix = "orca_corrected_snow_simple",                     label = "ORCA CADV",           years_from_end = 5),
    # (prefix = "orca_corrected_snow_cb0.12_ksymm500",            label = "ORCA Redi500",        years_from_end = 5),
    (prefix = "orca_corrected_snow_bih50days",           label = "ORCA CATKE",          years_from_end = 5),
    # (prefix = "orca_ncar_snow_cb0.15_bih50days",                label = "ORCA NCAR LowDiss",   years_from_end = 5),
    # (prefix = "halfdegree_corrected_snow_cb0.01_kskew0_ksymm0", label = "Half Degree CATKE",   years_from_end = 5),
    # (prefix = "orca_corrected_snow_cb0.06_kskew0_ksymm0",       label = "ORCA NOGM",           years_from_end = 5),
    (prefix = "orca_corrected_snow_kskew1000_ksymm1000_bih50days", label = "ORCA GM1000",         years_from_end = 5),
]

output_dir = length(ARGS) >= 1 ? ARGS[1] : "figures"

# ══════════════════════════════════════════════════════════════
# Infrastructure
# ══════════════════════════════════════════════════════════════

const HERE = @__DIR__
include(joinpath(HERE, "visualize", "common.jl"))
include(joinpath(HERE, "visualize", "cache.jl"))

# ══════════════════════════════════════════════════════════════
# Figure registry: (number, file basename, function symbol)
# ══════════════════════════════════════════════════════════════

const FIG_REGISTRY = [
    (n =  1, file = "fig01_sst_bias.jl",                  fn = :fig01),
    (n =  2, file = "fig02_sss_bias.jl",                  fn = :fig02),
    (n =  3, file = "fig03_ssh.jl",                       fn = :fig03),
    (n =  4, file = "fig04_mld.jl",                       fn = :fig04),
    (n =  5, file = "fig05_seaice_conc.jl",               fn = :fig05),
    (n =  6, file = "fig06_seaice_conc_bias.jl",          fn = :fig06),
    (n =  7, file = "fig07_surface_fluxes.jl",            fn = :fig07),
    (n =  8, file = "fig08_wind_stress.jl",               fn = :fig08),
    (n =  9, file = "fig09_ssh_variance.jl",              fn = :fig09),
    (n = 10, file = "fig10_sie.jl",                       fn = :fig10),
    (n = 11, file = "fig11_sia.jl",                       fn = :fig11),
    (n = 12, file = "fig12_arctic_volume.jl",             fn = :fig12),
    (n = 13, file = "fig13_sia_timeseries.jl",            fn = :fig13),
    (n = 14, file = "fig14_arctic_volume_timeseries.jl",  fn = :fig14),
    (n = 15, file = "fig15_ke.jl",                        fn = :fig15),
    (n = 16, file = "fig16_drift.jl",                     fn = :fig16),
    (n = 17, file = "fig17_profiles.jl",                  fn = :fig17),
    (n = 18, file = "fig18_zonal_mean.jl",                fn = :fig18),
    (n = 19, file = "fig19_zonal_drift.jl",               fn = :fig19),
    (n = 20, file = "fig20_mld_zonal_mean.jl",            fn = :fig20),
    (n = 21, file = "fig21_TS_drift_heatmap.jl",          fn = :fig21),
    (n = 22, file = "fig22_strait_transports.jl",         fn = :fig22),
    (n = 23, file = "fig23_amoc.jl",                      fn = :fig23),
    (n = 24, file = "fig24_near_surface_currents.jl",     fn = :fig24),
    (n = 25, file = "fig25_equatorial_undercurrent.jl",   fn = :fig25),
    (n = 26, file = "fig26_amoc_rapid.jl",                fn = :fig26),
]

# ══════════════════════════════════════════════════════════════
# Selection: parse FIG env var into a sorted, unique Vector{Int}
# Accepts "all" / unset → all figs; "4" → [4]; "1,4,7" → [1,4,7];
# "14-19" → [14..19]; mixed: "1,3,9-12,21".
# ══════════════════════════════════════════════════════════════

function parse_fig_selection(spec::AbstractString, all_ns::Vector{Int})
    s = strip(spec)
    (isempty(s) || lowercase(s) == "all") && return sort(unique(all_ns))
    out = Int[]
    for token in split(s, ',')
        t = strip(token)
        isempty(t) && continue
        if occursin('-', t)
            parts = split(t, '-')
            length(parts) == 2 || error("Invalid FIG range: '$t'")
            lo = parse(Int, strip(parts[1]))
            hi = parse(Int, strip(parts[2]))
            append!(out, lo:hi)
        else
            push!(out, parse(Int, t))
        end
    end
    return sort(unique(out))
end

# ══════════════════════════════════════════════════════════════
# Build per-case caches (cheap — no data loaded yet)
# Pre-include every fig file so `figNN` symbols are always defined.
# ══════════════════════════════════════════════════════════════

labels = [c.label for c in cases]
caches = Dict(c.label => CaseCache(c) for c in cases)

const FIGURES_DIR = joinpath(HERE, "visualize", "figures")

for entry in FIG_REGISTRY
    include(joinpath(FIGURES_DIR, entry.file))
end

# Convenience: render a single figure (or a list) by number from the REPL.
#     julia> render_figures(4)
#     julia> render_figures([1, 4, 7])
#     julia> render_figures([1, 4, 7]; theme = :dark)   # transparent bg, white axes
render_figures(n::Integer; kw...)                          = render_figures((n,); kw...)
render_figures(ns::AbstractVector{<:Integer}; kw...)       = render_figures(Tuple(ns); kw...)
render_figures(ns::AbstractRange{<:Integer}; kw...)        = render_figures(Tuple(ns); kw...)
function render_figures(ns::Tuple{Vararg{Integer}}; theme::Symbol = :light)
    with_theme(presentation_theme(theme)) do
        for entry in FIG_REGISTRY
            entry.n in ns || continue
            @info "Figure $(entry.n): $(entry.file)"
            getfield(@__MODULE__, entry.fn)(caches, labels, cases)
        end
    end
end

# ══════════════════════════════════════════════════════════════
# Auto-render only when invoked as a script (not from the REPL).
# ══════════════════════════════════════════════════════════════

if !isinteractive()
    selection = parse_fig_selection(get(ENV, "FIG", ""), [r.n for r in FIG_REGISTRY])
    theme = Symbol(lowercase(get(ENV, "THEME", "light")))
    @info "Rendering figures: $selection (theme = :$theme)"
    render_figures(selection; theme)
    @info "All requested figures saved to $output_dir"
else
    @info """
    visualize_omip.jl loaded.
      cases  = $(length(cases)) cases
      caches = pre-built (no data loaded yet)
    Render a figure with e.g.:  fig04(caches, labels, cases)
    Or render several with:     render_figures([1, 4, 17])
    """
end
