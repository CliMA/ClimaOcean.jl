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

"""
Perform `exp` only on rank 0, otherwise know as "root" rank.
Other ranks will wait for the root rank to finish before continuing
"""
macro root(exp)
    if MPI.Initialized() && MPI.Comm_rank(MPI.COMM_WORLD) != 0
        command = quote
          global_barrier()
        end
    else
        command = quote
           $exp
           global_barrier()
        end
    end

    return esc(command)
end

""" Distribute a `for` loop among ranks """
macro distribute(exp)
    if exp.head != :for
        error("The `@distribute` macro expects a `for` loop")
    end

    mpi_initialized = MPI.Initialized()
    if !mpi_initialized
        return esc(exp)
    end

    rank     = MPI.Comm_rank(MPI.COMM_WORLD)
    nprocs   = MPI.Comm_size(MPI.COMM_WORLD)
    iterable = exp.args[1].args[2]
    variable = exp.args[1].args[1]
    forbody  = exp.args[2]

    new_loop = quote
        for (counter, $variable) in enumerate($iterable)
            if $variable % $nprocs == $rank
                $forbody
            end
        end
        global_barrier()
    end

    return esc(new_loop)
end
