module DataWrangling

using Oceananigans
using Downloads
using Printf
using Downloads

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Grids: node
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: interpolate
using Oceananigans: pretty_filesize, location
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
            elapsed = 1e-9 * (time_ns() - download_start_time[])
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

function save_field_time_series!(fts; path, name, overwrite_existing=false)
    overwrite_existing && rm(path; force=true)

    times = on_architecture(CPU(), fts.times)
    grid  = on_architecture(CPU(), fts.grid)
    
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

include("metadata.jl")
include("inpaint_mask.jl")
include("JRA55/JRA55.jl")
include("ECCO/ECCO.jl")

using .ECCO
using .JRA55

end # module
