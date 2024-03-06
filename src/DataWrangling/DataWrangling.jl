module DataWrangling

using Oceananigans
using Downloads
using Printf

using Oceananigans.Architectures: architecture
using Oceananigans.Grids: node
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: interpolate
using Oceananigans: pretty_filesize
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

next_fraction = Ref(0.0)
download_start_time = Ref(time_ns())

"""
    download_progress(total, now; filename="")
"""
function download_progress(total, now; filename="")
    messages = 10

    if total > 0 
        fraction = now / total

        if fraction < 1 / messages && next_fraction[] == 0
            @info @sprintf("Downloading %s (size: %s)...", filename, pretty_filesize(total))
            next_fraction[] = 1 / messages
            download_start_time[] = time_ns()
        end

        if fraction > next_fraction[]
            elapsed = time_ns() - download_start_time[]
            msg = @sprintf(" ... downloaded %s (%d%% complete, %s)", pretty_filesize(now),
                           100fraction, prettytime(elapsed))
            @info msg
            next_fraction[] = next_fraction[] + 1 / messages
        end
    else
        if now > 0 && next_fraction[] == 0
            @info "Downloading $filename..."
            next_fraction[] = 1 / messages
            download_start_time[] = time_ns()
        end
    end

    return nothing
end

#####
##### FieldTimeSeries utilities
#####

# TODO: move this to Oceananigans

@kernel function _interpolate_field_time_series!(target_fts, target_grid, target_location,
                                                 source_fts, source_grid, source_location)

    # 4D index, cool!
    i, j, k, n = @index(Global, NTuple)

    source_field = view(source_fts, :, :, :, n)
    target_node = node(i, j, k, target_grid, target_location...)

    @inbounds target_fts[i, j, k, n] = interpolate(target_node, source_field, source_location, source_grid)
end

function interpolate_field_time_series!(target_fts, source_fts)
    @assert target_fts.times == source_fts.times
    times = target_fts.times
    Nt = length(times)

    target_grid = target_fts.grid
    source_grid = source_fts.grid

    @assert architecture(target_grid) == architecture(source_grid)
    arch = architecture(target_grid)

    # Make locations
    source_location = Tuple(L() for L in location(source_fts))
    target_location = Tuple(L() for L in location(target_fts))

    launch!(arch, target_grid, size(target_fts),
            _interpolate_field_time_series!,
            target_fts.data, target_grid, target_location,
            source_fts.data, source_grid, source_location)

    fill_halo_regions!(target_fts)

    return nothing
end

function save_field_time_series!(fts; path, name, overwrite_existing=false)
    overwrite_existing && rm(path; force=true)
    times = fts.times
    grid = fts.grid
    LX, LY, LZ = location(fts)
    ondisk_fts = FieldTimeSeries{LX, LY, LZ}(grid, times;
                                             backend = OnDisk(), path, name)

    Nt = length(times)
    for n = 1:Nt
        fill_halo_regions!(fts[n])
        set!(ondisk_fts, fts[n], n) 
    end

    return nothing
end

include("inpaint_mask.jl")
include("JRA55.jl")
include("ECCO2.jl")

using .JRA55
using .ECCO2

end # module
