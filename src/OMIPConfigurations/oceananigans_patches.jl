#####
##### Oceananigans monkey-patches and JLD2 split-file readers
#####
#
# Anything in this file is a deliberate override of upstream Oceananigans
# behaviour. Keep new patches here so they are easy to find, audit, and
# eventually upstream / remove. Each patch should carry a comment that says
# what it overrides, the upstream version it was written against, and why
# the override is needed.

using JLD2
using Oceananigans
using Oceananigans.Solvers
using Oceananigans.Solvers: ZDirection
using Oceananigans.Utils: worksize
using Oceananigans.Architectures: architecture
using Oceananigans.OutputReaders: SplitFilePath, InMemoryFTS, InMemory, Prefetched,
                                  time_indices, file_and_local_index
import Oceananigans.Fields: set!

#####
##### InMemory FieldTimeSeries split-file (`..._partN.jld2`) support
#####
#
# Oceananigans 0.107.x bug (`OutputReaders/field_time_series.jl`, ~line
# 924-930): the inner `FieldTimeSeries` constructor only builds a
# `SplitFilePath` when `backend isa OnDisk`. With an `InMemory` backend on
# split output (`..._part1.jld2`, `..._part2.jld2`, …), `fts.path` collapses
# to a single part file, producing "No data found for time …" warnings (and
# stale/zero data) on later snapshot reads.
#
# The patches below: (a) extend `set!` on an `InMemoryFTS` to iterate
# the per-part files referenced by a `SplitFilePath`, and (b) override the
# `FieldTimeSeries(path, name; …)` constructor so that an `InMemory`
# backend whose `fts.path` is a single file is detected as a split set and
# re-wrapped with the correct `SplitFilePath`.
#
# The helpers (`jld2_output_part_paths`, etc.) are factored out so this file
# and `scripts/visualize/common.jl` (and any future call site) can share a
# single source of truth for reading split Oceananigans JLD2 outputs.

"""
    jld2_output_part_paths(path)

Resolve `path` (Oceananigans-style JLD2 stem, with or without `.jld2`) to
the list of on-disk part files that hold its time series: either a
single-element vector `[abspath(path)]` if that file exists, or sorted
`…_partN.jld2` paths in `dirname(path)`. Returns an empty vector if nothing
matches (same convention as split detection for `FieldTimeSeries`).
"""
function jld2_output_part_paths(path::AbstractString)
    ap = abspath(path)
    isfile(ap) && return String[ap]
    base = endswith(ap, ".jld2") ? ap[1:end-5] : ap
    dir  = isempty(dirname(base)) ? "." : dirname(base)
    pat  = Regex("^" * Base.escape_string(basename(base)) * "_part(\\d+)\\.jld2\$")
    files = filter(f -> occursin(pat, f), readdir(dir))
    isempty(files) && return String[]
    sort!(files, by = f -> parse(Int, match(pat, f).captures[1]))
    return String[joinpath(dir, f) for f in files]
end

"""
    with_jld2(fn, path; reader_kw = NamedTuple())

Open `path` with JLD2, run `fn(file)`, guarantee close.
"""
function with_jld2(fn, path::AbstractString; reader_kw = NamedTuple())
    jf = JLD2.jldopen(path; reader_kw...)
    try
        return fn(jf)
    finally
        close(jf)
    end
end

"""
    jld2_parts(path)

Resolve `path` to its on-disk part files (single or split), erroring if
nothing matches. Used by every metadata reader below.
"""
function jld2_parts(path::AbstractString)
    parts = jld2_output_part_paths(path)
    isempty(parts) && error("No JLD2 output at path '$path' (single file or split parts).")
    return parts
end

# Per-part memo: `fn(file)` runs at most once per absolute part path, its
# return value cached in `memo`. Each derived diagnostic that validates a
# disk cache hits this memo and is therefore O(1) after the first read.
const JLD2_NT_PER_PART    = Dict{String, Int}()
const JLD2_TIMES_PER_PART = Dict{String, Vector{Float64}}()

memoize_jld2_part(fn, memo, path::AbstractString; reader_kw = NamedTuple()) =
    get!(memo, path) do
        with_jld2(fn, path; reader_kw)
    end

"""
    total_jld2_timeseries_snapshot_count(path; reader_kw = NamedTuple())

Total number of time indices for an Oceananigans JLD2 output stem `path`
(single file or split `_partN.jld2` parts), summed across parts. Reads
JLD2 metadata only — no `FieldTimeSeries` is built. Per-part counts are
memoized for the session.
"""
function total_jld2_timeseries_snapshot_count(path::AbstractString; reader_kw = NamedTuple())
    n = 0
    for p in jld2_parts(path)
        n += memoize_jld2_part(JLD2_NT_PER_PART, p; reader_kw) do jf
            length(keys(jf["timeseries/t"]))
        end
    end
    return n
end

"""
    total_jld2_timeseries_times(path; reader_kw = NamedTuple())

Concatenated time-coordinate vector for an Oceananigans JLD2 output stem
`path`. Reads `timeseries/t/<iter>` per part, sorted by iteration, and
concatenates across parts. No `FieldTimeSeries` construction. Per-part
time vectors are memoized for the session.
"""
function total_jld2_timeseries_times(path::AbstractString; reader_kw = NamedTuple())
    times = Float64[]
    for p in jld2_parts(path)
        ts = memoize_jld2_part(JLD2_TIMES_PER_PART, p; reader_kw) do jf
            iterations = sort!(parse.(Int, collect(keys(jf["timeseries/t"]))))
            return Float64[jf["timeseries/t/$it"] for it in iterations]
        end
        append!(times, ts)
    end
    return times
end

"""
    last_jld2_timeseries_time(path; reader_kw = NamedTuple())

Latest snapshot time for an Oceananigans JLD2 output stem `path`. Opens
only the last part file and reads only the highest-iteration entry of
`timeseries/t`. Avoids the per-iteration JLD2 lookup loop in
`total_jld2_timeseries_times`, which is O(N) and slow on long runs over
network filesystems.
"""
function last_jld2_timeseries_time(path::AbstractString; reader_kw = NamedTuple())
    last_part = last(jld2_parts(path))
    return with_jld2(last_part; reader_kw) do jf
        ks = keys(jf["timeseries/t"])
        max_iter = maximum(parse(Int, k) for k in ks if !isnothing(tryparse(Int, k)))
        return Float64(jf["timeseries/t/$max_iter"])
    end
end

"""
    total_jld2_serialized_grid(path, name; reader_kw = NamedTuple())

Read the grid serialized inside an Oceananigans JLD2 output stem `path`
for variable `name`, without instantiating a `FieldTimeSeries`. Defers to
`Oceananigans.OutputReaders.load_serialized_grid`, which handles
single-grid (`serialized/grid`) and multi-grid output. The grid is
invariant across split parts, so only the first part is opened.
"""
total_jld2_serialized_grid(path::AbstractString, name::AbstractString; reader_kw = NamedTuple()) =
    with_jld2(first(jld2_parts(path)); reader_kw) do jf
        Oceananigans.OutputReaders.load_serialized_grid(jf, name)
    end

"""
    total_jld2_scalar_timeseries(path, name; reader_kw = NamedTuple())

Read a scalar (1×1×1) variable's full time series for the JLD2 output stem
`path` directly from `timeseries/<name>/<iter>` metadata, sorted by
iteration and concatenated across split parts. Avoids building a
`FieldTimeSeries` and the per-snapshot buffer slides that imposes — ~100×
faster than `[interior(fts[n])[1] for n in ...]` for long time series of
`tosga`/`soga`-style global means.
"""
function total_jld2_scalar_timeseries(path::AbstractString, name::AbstractString;
                                      reader_kw = NamedTuple())
    values = Float64[]
    for p in jld2_parts(path)
        chunk = with_jld2(p; reader_kw) do jf
            ks = collect(keys(jf["timeseries/$name"]))
            iterations = sort!([parse(Int, k) for k in ks if !isnothing(tryparse(Int, k))])
            return Float64[jf["timeseries/$name/$it"][1] for it in iterations]
        end
        append!(values, chunk)
    end
    return values
end

function detect_split_file_path(path::AbstractString, reader_kw)
    isfile(path) && return nothing
    parts = jld2_output_part_paths(path)
    isempty(parts) && return nothing
    nper = Int[total_jld2_timeseries_snapshot_count(p; reader_kw) for p in parts]
    return SplitFilePath(parts, cumsum(nper))
end

location_types(::Oceananigans.OutputReaders.FieldTimeSeries{LX, LY, LZ}) where {LX, LY, LZ} = (LX, LY, LZ)

rebuild_backend_with_path(backend, new_path) = backend

function rebuild_backend_with_path(backend::Prefetched, new_path)
    old_buf = getfield(backend, :buffer_fts)
    BLX, BLY, BLZ = location_types(old_buf)
    new_buf = Oceananigans.OutputReaders.FieldTimeSeries{BLX, BLY, BLZ}(
        old_buf.data, old_buf.grid, old_buf.backend, old_buf.boundary_conditions,
        old_buf.indices, old_buf.times, new_path, old_buf.name,
        old_buf.time_indexing, old_buf.reader_kw)
    return Prefetched(backend.base_backend, backend.pending, new_buf, backend.next_start)
end

function rebuild_fts_with_path(fts, new_path)
    LX, LY, LZ = location_types(fts)
    new_backend = rebuild_backend_with_path(fts.backend, new_path)
    return Oceananigans.OutputReaders.FieldTimeSeries{LX, LY, LZ}(
        fts.data, fts.grid, new_backend, fts.boundary_conditions,
        fts.indices, fts.times, new_path, fts.name,
        fts.time_indexing, fts.reader_kw)
end

# Patch 1: iterate per-part files when an InMemory FTS is `set!` from a
# `SplitFilePath` (Oceananigans does this only for `OnDisk` in 0.107.x).
function set!(fts::InMemoryFTS, sfp::SplitFilePath, name::String = fts.name;
              warn_missing_data = false, kwargs...)
    idxs = time_indices(fts)
    Ntot = last(sfp.cumulative_length)
    needed = String[]
    for n in idxs
        (n < 1 || n > Ntot) && continue
        file_path, _ = file_and_local_index(sfp, n)
        file_path ∉ needed && push!(needed, file_path)
    end
    for p in needed
        set!(fts, p, name; warn_missing_data, kwargs...)
    end
    return nothing
end

# Patch 2: detect split sets when the user passes a single stem path with an
# `InMemory` backend, and rewrap the FTS so its `path` is a `SplitFilePath`.
function Oceananigans.OutputReaders.FieldTimeSeries(path::String, name::String;
                                                    backend = InMemory(),
                                                    reader_kw = NamedTuple(),
                                                    kwargs...)
    fts = invoke(Oceananigans.OutputReaders.FieldTimeSeries,
                 Tuple{String, Vararg{Any}},
                 path, name; backend, reader_kw, kwargs...)
    if backend isa InMemory && !(fts.path isa SplitFilePath)
        sfp = detect_split_file_path(path, reader_kw)
        sfp === nothing || (fts = rebuild_fts_with_path(fts, sfp))
    end
    return fts
end
