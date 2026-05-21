using Dates
using Oceananigans.DistributedComputations: @root,
                                            mpi_initialized,
                                            mpi_rank,
                                            global_communicator

# JRA55 shortnames used in filenames (subset actually loaded by OMIP)
const JRA55_SHORTNAMES = ["tas", "huss", "psl", "uas", "vas", "rlds", "rsds", "prra", "prsn", "friver", "licalvf"]

"""
    setup_staging_directory(source_dir, staging_dir)

Populate `staging_dir` with symlinks to every `.nc` file in `source_dir`.
Reads go through symlinks (slow) until files are staged with `stage_jra55_year!`.
Also:
  - sweeps any leftover `*.tmp` from a killed prior run (half-written copies
    that would otherwise be mistaken for real copies);
  - replaces any real `.nc` file whose size doesn't match the source with
    a fresh symlink. A size mismatch means a prior run crashed mid-`cp` and
    left a truncated file — before this heal, HDF5 reads would fail with
    `truncated file: eof = X, stored_eof = Y`, and `stage_jra55_year!`
    would skip it (it only replaces symlinks).
"""
function setup_staging_directory(source_dir, staging_dir)
    # Filesystem mutations only on rank 0; other ranks barrier at the end of
    # @root, so by the time they return all symlinks/heals are visible.
    @root begin
        mkpath(staging_dir)
        for leftover in filter(f -> endswith(f, ".nc.tmp"), readdir(staging_dir; join=true))
            rm(leftover; force=true)
        end
        for src in filter(f -> endswith(f, ".nc"), readdir(source_dir; join=true))
            dst = joinpath(staging_dir, basename(src))
            if islink(dst)
                continue                                         # healthy symlink, leave alone
            elseif isfile(dst)                                    # real file — validate size
                if filesize(dst) != filesize(src)
                    @warn "setup_staging_directory: size mismatch at $dst ($(filesize(dst)) vs source $(filesize(src))); replacing with symlink"
                    rm(dst; force=true)
                    symlink(src, dst)
                end
            elseif !ispath(dst)                                   # nothing there (incl. broken symlink handled by islink above)
                symlink(src, dst)
            end
        end
    end
    return staging_dir
end

# Replace `dst` atomically with whatever `make_new!(tmp)` produces at `tmp`.
# `rename(2)` is atomic on the same filesystem, so concurrent readers either
# see the old `dst` (symlink or previous real copy) or the new one — never a
# half-written file. Readers holding an fd to the old inode keep reading it
# correctly; the kernel keeps the inode alive until they close.
function atomic_replace!(dst, make_new!)
    tmp = dst * ".tmp"
    isfile(tmp) && rm(tmp; force=true)     # stale tmp from a crash — drop it
    make_new!(tmp)
    mv(tmp, dst; force=true)
    return dst
end

# Pure file-system work. No MPI calls — safe to run from a background thread
# launched on rank 0. Concurrency with the main task is mediated by `atomic_replace!`
# (rename(2) is atomic on the same filesystem), so readers see either the old
# symlink or the new real copy, never a half-written file.
function _stage_jra55_year_files!(source_dir, staging_dir, year)
    year_str = string(year)
    for name in JRA55_SHORTNAMES
        for dst in filter(f -> contains(f, name) && contains(f, year_str) && endswith(f, ".nc"), readdir(staging_dir; join=true))
            if islink(dst)
                src = joinpath(source_dir, basename(dst))
                atomic_replace!(dst, tmp -> cp(src, tmp))
                @debug "Staged $(basename(dst)) to scratch"
            end
        end
    end
    return nothing
end

function _unstage_jra55_year_files!(source_dir, staging_dir, year)
    year_str = string(year)
    for name in JRA55_SHORTNAMES
        for dst in filter(f -> contains(f, name) && contains(f, year_str) && endswith(f, ".nc"), readdir(staging_dir; join=true))
            if isfile(dst) && !islink(dst)
                src = joinpath(source_dir, basename(dst))
                atomic_replace!(dst, tmp -> symlink(src, tmp))
                @debug "Unstaged $(basename(dst)) from scratch"
            end
        end
    end
    return nothing
end

"""
    stage_jra55_year!(source_dir, staging_dir, year)

Replace symlinks in `staging_dir` with real copies for all JRA55 files
matching `year`. Skips files that are already real copies. The replacement
is atomic — a partial copy is never visible at `dst`, so concurrent
`PrefetchingBackend` readers on background threads cannot race with the `cp`.
"""
function stage_jra55_year!(source_dir, staging_dir, year)
    @root _stage_jra55_year_files!(source_dir, staging_dir, year)
    return nothing
end

"""
    unstage_jra55_year!(source_dir, staging_dir, year)

Remove real copies for `year` from `staging_dir` and restore symlinks
to `source_dir`. Uses the same atomic swap as staging so a concurrent
reader never sees a missing file at the swap point.
"""
function unstage_jra55_year!(source_dir, staging_dir, year)
    @root _unstage_jra55_year_files!(source_dir, staging_dir, year)
    return nothing
end

"""
    JRA55DataStagingCallback(; source_dir, staging_dir, start_date, async = true)

Return a simulation callback that dynamically stages JRA55 yearly files
from `source_dir` (slow disk) to `staging_dir` (fast scratch).

At each invocation the callback:
  1. Reaps any background staging tasks that have finished
  2. If the *current* year's data isn't ready yet, blocks until it is
     (safety fallback — should be rare if the callback fires often enough)
  3. Spawns background copies for the current and next year's files
     (if they aren't already on scratch or in flight)
  4. Removes files from two or more years ago to free space

Each year of JRA55 data is ~15–25 GB, so scratch holds at most ~50 GB at any time.

# Async vs sync

With `async = true` (default), the actual `cp` is launched on a Julia
background thread via `Threads.@spawn` while the main task returns
immediately. With Δt = 30 min, ~58 min of compute fits between year
transitions, so a 20-min copy overlaps fully with compute. **Start Julia
with `JULIA_NUM_THREADS ≥ 2`** for this to actually parallelise — with one
thread, `@spawn` just runs the copy synchronously.

The atomic `mv(tmp, dst)` inside `atomic_replace!` is the only point where
concurrent readers can observe the file change. They see either the old
symlink (slow read) or the new real copy (fast read), never a partial
file. On MPI runs, the cp itself runs only on rank 0 in a background
thread; cross-rank consistency comes from filesystem atomicity rather than
an explicit barrier, so no MPI calls happen on the background thread.

Set `async = false` to fall back to the original blocking behaviour.
"""
function JRA55DataStagingCallback(; source_dir, staging_dir, start_date, async = true)

    # "We've already issued a stage request for this year" — kept identical
    # across ranks so every `@root` / `@root @info` below fires symmetrically
    # (the macro contains an `MPI.Barrier`, so asymmetric entry deadlocks).
    requested_years = Set{Int}()

    # In-flight per-year background tasks. Only rank 0 ever has entries.
    staging_tasks   = Dict{Int, Task}()

    # Reap background tasks that have finished. No MPI calls, no `@root`:
    # other ranks see an empty dict and silently no-op.
    function reap_completed_tasks!()
        for (y, task) in collect(staging_tasks)
            if istaskdone(task)
                if istaskfailed(task)
                    err = try
                        Base.task_result(task)
                    catch
                        nothing
                    end
                    @warn "Async staging for year $y failed; reads will fall back to the source symlink." exception = err
                end
                delete!(staging_tasks, y)
            end
        end
    end

    # If a task for `y` is still running and the simulation needs the file
    # now, block on it. Only rank 0 ever has tasks, so this branch is
    # rank-0-only; we deliberately use plain `@info` (not `@root @info`) so
    # no barrier fires. Other ranks continue and will re-synchronise with
    # rank 0 at the next halo-exchange MPI call.
    function ensure_year_ready!(y)
        if haskey(staging_tasks, y)
            @info "Background staging of year $y not done yet; rank 0 is blocking on the copy."
            try
                wait(staging_tasks[y])
            catch e
                @warn "Async staging of year $y errored mid-wait" exception = e
            end
            delete!(staging_tasks, y)
        end
    end

    # On rank 0, kick off the cp on a background thread; on every other
    # rank, do nothing — the shared filesystem (with `atomic_replace!`)
    # makes the new file visible without an explicit MPI sync.
    function spawn_stage_year!(y)
        if !async
            stage_jra55_year!(source_dir, staging_dir, y)
            return
        end
        if !mpi_initialized() || mpi_rank(global_communicator()) == 0
            staging_tasks[y] = Threads.@spawn _stage_jra55_year_files!(source_dir, staging_dir, y)
        end
    end

    function stage_forcing_data!(sim)
        current_year = year(start_date + Second(round(Int, sim.model.clock.time)))
        needed_years = (current_year, current_year + 1)

        reap_completed_tasks!()

        # Safety fallback: if the current simulated year hasn't finished
        # copying yet, rank 0 blocks until it has. (Should be rare — the
        # 20 min copy fits inside ~58 min of compute per simulated year.)
        ensure_year_ready!(current_year)

        # Schedule staging for the years we need. `requested_years` is
        # updated *before* the spawn on every rank, so the `@root @info`
        # below fires identically on rank 0 and on every other rank.
        for y in needed_years
            if y ∉ requested_years
                push!(requested_years, y)
                @root @info "$(async ? "Spawning async" : "Staging") JRA55 data for year $y → $staging_dir"
                spawn_stage_year!(y)
            end
        end

        # Unstage years we no longer need. `requested_years` is consistent
        # across ranks, so the `@root` inside `unstage_jra55_year!` fires
        # on the same set of years on every rank.
        for y in collect(requested_years)
            if y < current_year - 1
                @root @info "Unstaging JRA55 data for year $y from $staging_dir"
                unstage_jra55_year!(source_dir, staging_dir, y)
                delete!(requested_years, y)
            end
        end
    end

    return stage_forcing_data!
end
