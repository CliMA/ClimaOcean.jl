using MPI

#####
##### Handle commands, typically downloading files
##### which should be executed by only one rank or distributed among ranks
#####


function global_barrier()
    if MPI.Initialized()
        MPI.Barrier(MPI.COMM_WORLD)
    end
end

function blocking_run(cmd)
    if MPI.Initialized() && MPI.Comm_rank(MPI.COMM_WORLD) != 0
        nothing # Do nothing!
    else
        run(cmd)
    end

    global_barrier()
    
    return nothing
end

function blocking_download(url, filepath; overwrite_existing = false, kw...)
    if overwrite_existing
        blocking_run(`rm $filepath`)
    end

    if !isfile(filepath) 
        if MPI.Initialized() && MPI.Comm_rank(MPI.COMM_WORLD) != 0
            nothing # Do nothing!
        else
            Downloads.download(url, filepath; kw...)
        end    
    end

    global_barrier()

    return nothing
end

# If only one command is provided, run it with rank 0
distributed_run(cmd) = blocking_run(cmd)

# If multiple commands are provided, distribute them among the ranks
function distributed_run(cmds::Union{AbstractVector, Tuple})
    if length(cmds) == 1 # just one command
        distributed_run(cmds[1])
    else # split among different ranks
        rank   = MPI.Comm_rank(MPI.COMM_WORLD)
        nprocs = MPI.Comm_size(MPI.COMM_WORLD)
        for (i, cmd) in enumerate(cmds)
            if i % nprocs == rank
                blocking_run(cmd)
            end
        end
    end
    global_barrier()

    return nothing
end
