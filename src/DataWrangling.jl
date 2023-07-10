module DataWrangling

using Downloads
using Printf

next_fraction = Ref(0.0)
download_start_time = Ref(time_ns())

function download_progress(total, now)
    messages = 10

    if total > 0 
        fraction = now / total

        if fraction < 1 / messages && next_fraction[] == 0
            @info @sprintf("Downloading bathymetry (file size: %s)...", pretty_filesize(total))
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
            @info "Downloading bathymetry..."
            next_fraction[] = 1 / messages
            download_start_time[] = time_ns()
        end
    end

    return nothing
end

end # module
