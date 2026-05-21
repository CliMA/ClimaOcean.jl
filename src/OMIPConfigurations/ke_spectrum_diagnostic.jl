using FFTW
using JLD2
using CUDA
using LinearAlgebra: mul!
using AbstractFFTs: plan_fft
using Oceananigans.Architectures: architecture, on_architecture, child_architecture, AbstractArchitecture, CPU
using Oceananigans.Grids: λnode, φnode, znode, Center
using Oceananigans.Fields: interior
using Oceananigans.Utils: prettytime, TimeInterval
using Oceananigans.Simulations: add_callback!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.DistributedComputations: Distributed, concatenate_local_sizes,
                                              reconstruct_global_grid
import MPI

# A 2-D rectangular patch in (i, j) at a list of vertical k-levels.
struct KEPatch
    name      :: String
    i_range   :: UnitRange{Int}
    j_range   :: UnitRange{Int}
    k_levels  :: Vector{Int}
    depths    :: Vector{Float64}
    λ_center  :: Float64
    φ_center  :: Float64
end

# NOTE: For distributed runs, `size(grid)` returns LOCAL sizes and λnode/φnode/znode
# are evaluated on the local subgrid. All index computation here MUST run on the
# reconstructed *global* CPU grid — see `add_ke_spectrum_diagnostic!`.

# Brute-force search over the global grid for the (i, j) cell whose (λ, φ) is closest
# to (λ_target, φ_target). Runs once at setup; cost ~1s on 3600×1800.
function nearest_grid_index(global_cpu_grid, λ_target, φ_target)
    Nx, Ny, Nz = size(global_cpu_grid)
    best_i, best_j, best_d = 0, 0, Inf
    for j in 1:Ny, i in 1:Nx
        λ = λnode(i, j, Nz, global_cpu_grid, Center(), Center(), Center())
        φ = φnode(i, j, Nz, global_cpu_grid, Center(), Center(), Center())
        dλ = mod(λ - λ_target + 180, 360) - 180
        dφ = φ - φ_target
        d = dλ * dλ + dφ * dφ
        if d < best_d
            best_d = d
            best_i, best_j = i, j
        end
    end
    return best_i, best_j
end

function nearest_k_indices(global_cpu_grid, target_depths)
    Nz = size(global_cpu_grid, 3)
    k_levels = Int[]
    actual_depths = Float64[]
    for d_target in target_depths
        best_k, best_err = 0, Inf
        for k in 1:Nz
            z = znode(1, 1, k, global_cpu_grid, Center(), Center(), Center())
            err = abs(-z - d_target)
            if err < best_err
                best_err = err
                best_k = k
            end
        end
        push!(k_levels, best_k)
        push!(actual_depths, -znode(1, 1, best_k, global_cpu_grid, Center(), Center(), Center()))
    end
    return k_levels, actual_depths
end

function build_patch(name, λ_center, φ_center, Nx_patch, Ny_patch, target_depths, global_cpu_grid)
    ic, jc = nearest_grid_index(global_cpu_grid, λ_center, φ_center)
    Δi, Δj = Nx_patch ÷ 2, Ny_patch ÷ 2
    i_range = (ic - Δi):(ic + Δi - 1)
    j_range = (jc - Δj):(jc + Δj - 1)
    Nx, Ny = size(global_cpu_grid, 1), size(global_cpu_grid, 2)
    (first(i_range) >= 1 && last(i_range) <= Nx) ||
        error("KE-spectrum patch $name: i_range = $i_range out of bounds (1..$Nx)")
    (first(j_range) >= 1 && last(j_range) <= Ny) ||
        error("KE-spectrum patch $name: j_range = $j_range out of bounds (1..$Ny)")
    k_levels, depths = nearest_k_indices(global_cpu_grid, target_depths)
    return KEPatch(name, i_range, j_range, k_levels, depths, λ_center, φ_center)
end

# Separable 2-D Hann window.
function hann_window(::Type{T}, Nx, Ny) where T
    wx = T[T(0.5) * (1 - cos(2 * T(π) * (i - 1) / (Nx - 1))) for i in 1:Nx]
    wy = T[T(0.5) * (1 - cos(2 * T(π) * (j - 1) / (Ny - 1))) for j in 1:Ny]
    return wx .* wy'
end

# Global (i, j) range owned by this rank, given pre-computed per-rank size vectors
# (computed once at setup so we don't repeat MPI.Allreduce calls per gather).
function rank_global_ranges(arch::Distributed, nx_per_rank, ny_per_rank)
    ri, rj, _ = arch.local_index
    i_glob = (1 + sum(nx_per_rank[1:ri-1])) : sum(nx_per_rank[1:ri])
    j_glob = (1 + sum(ny_per_rank[1:rj-1])) : sum(ny_per_rank[1:rj])
    return i_glob, j_glob
end

# Each rank fills its intersection of `patch` at level `k` into `host_buf`; MPI.Allreduce!
# sums across ranks so every rank holds the full (Nx_patch, Ny_patch) patch on the host.
function gather_patch_to_host!(host_buf::AbstractMatrix, field, patch::KEPatch, k::Int,
                                arch::Distributed, nx_per_rank, ny_per_rank)
    fill!(host_buf, 0)
    i_glob, j_glob = rank_global_ranges(arch, nx_per_rank, ny_per_rank)

    i_lo = max(first(patch.i_range), first(i_glob))
    i_hi = min(last(patch.i_range),  last(i_glob))
    j_lo = max(first(patch.j_range), first(j_glob))
    j_hi = min(last(patch.j_range),  last(j_glob))

    if i_lo <= i_hi && j_lo <= j_hi
        i_loc = (i_lo - first(i_glob) + 1):(i_hi - first(i_glob) + 1)
        j_loc = (j_lo - first(j_glob) + 1):(j_hi - first(j_glob) + 1)
        i_buf = (i_lo - first(patch.i_range) + 1):(i_hi - first(patch.i_range) + 1)
        j_buf = (j_lo - first(patch.j_range) + 1):(j_hi - first(patch.j_range) + 1)
        slice = Array(view(interior(field), i_loc, j_loc, k))
        @views host_buf[i_buf, j_buf] .= slice
    end

    MPI.Allreduce!(host_buf, +, arch.communicator)
    return host_buf
end

# Single-rank fallback (CPU or single-GPU): direct device-to-host copy of the patch.
function gather_patch_to_host!(host_buf::AbstractMatrix, field, patch::KEPatch, k::Int,
                                ::AbstractArchitecture, _, _)
    slice = Array(view(interior(field), patch.i_range, patch.j_range, k))
    copyto!(host_buf, slice)
    return host_buf
end

mutable struct KESpectrumDiagnostic{Fu, Fv, P, R, C, B, W}
    u_field        :: Fu
    v_field        :: Fv
    patches        :: Vector{KEPatch}
    Nx_patch       :: Int
    Ny_patch       :: Int
    window         :: W                  # device-resident, real
    fft_plan       :: P                  # device-side, complex→complex
    host_buf_u     :: Matrix{Float64}    # host scratch for gather (real)
    host_buf_v     :: Matrix{Float64}
    dev_real_u     :: R                  # device real, post host→device copy
    dev_real_v     :: R
    dev_cplx_u     :: C                  # device complex (windowed FFT input)
    dev_cplx_v     :: C
    dev_fft_u      :: C                  # device complex (FFT output)
    dev_fft_v      :: C
    psd_sum        :: Vector{Vector{B}}  # per (patch, level) device real accumulator
    counters       :: Vector{Int}
    nx_per_rank    :: Vector{Int}        # cached at construction (Allreduce-free per call)
    ny_per_rank    :: Vector{Int}
    output_path    :: String
    flush_interval :: Float64
    start_time     :: Float64            # don't accumulate until sim.clock.time >= start_time
    last_flush     :: Base.RefValue{Float64}
    is_root        :: Bool
    started        :: Base.RefValue{Bool}
end

function build_ke_spectrum_diagnostic(model, patches, Nx_patch, Ny_patch,
                                       output_path, flush_interval, start_time)
    u = model.ocean.model.velocities.u
    v = model.ocean.model.velocities.v
    grid = model.ocean.model.grid
    arch = architecture(grid)
    child = arch isa Distributed ? child_architecture(arch) : arch
    is_root = arch isa Distributed ? all(arch.local_index .== 1) : true

    # Cache per-rank size tuples once. `concatenate_local_sizes` does an MPI.Allreduce
    # internally; calling it per gather would mean tens of redundant Allreduces per cycle.
    if arch isa Distributed
        sizes = concatenate_local_sizes(size(grid), arch)
        nx_per_rank = collect(sizes[1])
        ny_per_rank = collect(sizes[2])
    else
        nx_per_rank = Int[size(grid, 1)]
        ny_per_rank = Int[size(grid, 2)]
    end

    window_dev = on_architecture(child, hann_window(Float64, Nx_patch, Ny_patch))

    dev_real_u = on_architecture(child, zeros(Float64,    Nx_patch, Ny_patch))
    dev_real_v = on_architecture(child, zeros(Float64,    Nx_patch, Ny_patch))
    dev_cplx_u = on_architecture(child, zeros(ComplexF64, Nx_patch, Ny_patch))
    dev_cplx_v = on_architecture(child, zeros(ComplexF64, Nx_patch, Ny_patch))
    dev_fft_u  = on_architecture(child, zeros(ComplexF64, Nx_patch, Ny_patch))
    dev_fft_v  = on_architecture(child, zeros(ComplexF64, Nx_patch, Ny_patch))

    fft_plan = plan_fft(dev_cplx_u)

    psd_sum = [
        [on_architecture(child, zeros(Float64, Nx_patch, Ny_patch)) for _ in patch.k_levels]
        for patch in patches
    ]
    counters = zeros(Int, length(patches))

    host_buf_u = zeros(Float64, Nx_patch, Ny_patch)
    host_buf_v = zeros(Float64, Nx_patch, Ny_patch)

    return KESpectrumDiagnostic(u, v, patches, Nx_patch, Ny_patch,
                                window_dev, fft_plan,
                                host_buf_u, host_buf_v,
                                dev_real_u, dev_real_v,
                                dev_cplx_u, dev_cplx_v,
                                dev_fft_u, dev_fft_v,
                                psd_sum, counters,
                                nx_per_rank, ny_per_rank,
                                output_path,
                                Float64(flush_interval), Float64(start_time),
                                Ref(0.0), is_root, Ref(false))
end

function (D::KESpectrumDiagnostic)(sim)
    sim_time = Float64(sim.model.clock.time)
    sim_time < D.start_time && return nothing
    if !D.started[]
        D.last_flush[] = sim_time
        D.started[] = true
    end

    grid = sim.model.ocean.model.grid
    arch = architecture(grid)

    for (p, patch) in enumerate(D.patches)
        for (l, k) in enumerate(patch.k_levels)
            gather_patch_to_host!(D.host_buf_u, D.u_field, patch, k, arch,
                                   D.nx_per_rank, D.ny_per_rank)
            gather_patch_to_host!(D.host_buf_v, D.v_field, patch, k, arch,
                                   D.nx_per_rank, D.ny_per_rank)
            copyto!(D.dev_real_u, D.host_buf_u)
            copyto!(D.dev_real_v, D.host_buf_v)
            @. D.dev_cplx_u = ComplexF64(D.dev_real_u * D.window)
            @. D.dev_cplx_v = ComplexF64(D.dev_real_v * D.window)
            mul!(D.dev_fft_u, D.fft_plan, D.dev_cplx_u)
            mul!(D.dev_fft_v, D.fft_plan, D.dev_cplx_v)
            @. D.psd_sum[p][l] += abs2(D.dev_fft_u) + abs2(D.dev_fft_v)
        end
        D.counters[p] += 1
    end

    if sim_time - D.last_flush[] >= D.flush_interval
        D.is_root && flush_spectrum!(D, sim_time)
        for arrs in D.psd_sum, m in arrs
            fill!(m, 0)
        end
        D.counters .= 0
        D.last_flush[] = sim_time
    end

    return nothing
end

function flush_spectrum!(D::KESpectrumDiagnostic, sim_time::Float64)
    window_norm_sq = sum(abs2, Array(D.window))
    JLD2.jldopen(D.output_path, "a+") do file
        for (p, patch) in enumerate(D.patches)
            D.counters[p] == 0 && continue
            for (l, _) in enumerate(patch.k_levels)
                psd_host = Array(D.psd_sum[p][l]) ./ (D.counters[p] * window_norm_sq)
                key = "$(patch.name)/level$(l)/t_$(sim_time)"
                file[key] = psd_host
            end
            meta_key = "$(patch.name)/_metadata"
            if !haskey(file, meta_key)
                file[meta_key] = (i_range = collect(patch.i_range),
                                  j_range = collect(patch.j_range),
                                  k_levels = patch.k_levels,
                                  depths = patch.depths,
                                  λ_center = patch.λ_center,
                                  φ_center = patch.φ_center)
            end
        end
        if !haskey(file, "_window/Nx")
            file["_window/Nx"] = D.Nx_patch
            file["_window/Ny"] = D.Ny_patch
            file["_window/type"] = "Hann"
            file["_window/norm_sq"] = window_norm_sq
        end
    end
    return nothing
end

# Default no-op for non-tenthdegree configurations.
add_ke_spectrum_diagnostic!(simulation, ::Val; kwargs...) = nothing

function add_ke_spectrum_diagnostic!(simulation, ::Val{:tenthdegree};
                                      output_dir = ".",
                                      filename_prefix = "tenthdegree",
                                      fft_interval = 1hours,
                                      flush_interval = 15days,
                                      start_time = 2 * 365days,
                                      Nx_patch = 256,
                                      Ny_patch = 256,
                                      target_depths = (0.0, 100.0, 500.0, 1000.0))
    grid = simulation.model.ocean.model.grid

    # Patch indices live in the *global* index space. `reconstruct_global_grid` is a
    # no-op for non-distributed grids, so this works for all configs.
    #
    # We unwrap to the underlying (non-immersed) grid first because Oceananigans'
    # `reconstruct_global_grid(::ImmersedBoundaryGrid)` has a bug at
    # `distributed_immersed_boundaries.jl:28`: it passes `size(grid)` (i.e.
    # `(Nx, Ny, Nz)`) to `construct_global_array` for the bottom_height gather,
    # which makes the global bottom_height a 3D `(Nx, Ny, Nz)` array. Then
    # `materialize_immersed_boundary` aborts with a BoundsError trying to assign
    # that 3D array into the 2D `Field{Center, Center, Nothing}` bottom_height.
    # `build_patch` only reads coordinates, never the immersed boundary, so we
    # don't need the IBG layer here at all.
    underlying = grid isa ImmersedBoundaryGrid ? grid.underlying_grid : grid
    global_grid     = reconstruct_global_grid(underlying)
    global_cpu_grid = on_architecture(CPU(), global_grid)

    patch_specs = (
        (name = "gulf_stream",        λ = 295.0, φ =  37.5),
        (name = "kuroshio",           λ = 160.0, φ =  37.5),
        (name = "acc_drake",          λ = 310.0, φ = -50.0),
        (name = "equatorial_pacific", λ = 210.0, φ =   0.0),
    )

    patches = [build_patch(s.name, s.λ, s.φ, Nx_patch, Ny_patch,
                            collect(target_depths), global_cpu_grid) for s in patch_specs]

    output_path = joinpath(output_dir, filename_prefix * "_ke_spectrum.jld2")
    diag = build_ke_spectrum_diagnostic(simulation.model, patches, Nx_patch, Ny_patch,
                                         output_path, flush_interval, start_time)

    add_callback!(simulation, diag, TimeInterval(fft_interval))

    diag.is_root && @info "KE spectrum diagnostic attached: " *
                          "$(length(patches)) patches × $(length(target_depths)) levels, " *
                          "starts at $(prettytime(start_time)), " *
                          "FFT every $(prettytime(fft_interval)), " *
                          "flush every $(prettytime(flush_interval))"

    return diag
end

add_ke_spectrum_diagnostic!(simulation, config::Symbol; kwargs...) =
    add_ke_spectrum_diagnostic!(simulation, Val(config); kwargs...)
