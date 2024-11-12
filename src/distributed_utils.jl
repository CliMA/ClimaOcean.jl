using MPI

#####
##### Handle commands, typically downloading files
##### which should be executed by only one rank or distributed among ranks
#####

# Utilities to make the macro work importing only ClimaOcean and not MPI
mpi_initialized() = MPI.Initialized()
mpi_rank()        = MPI.Comm_rank(MPI.COMM_WORLD)
mpi_size()        = MPI.Comm_size(MPI.COMM_WORLD)
global_barrier()  = mpi_initialized() ? MPI.Barrier(MPI.COMM_WORLD) : nothing

"""
    @root exs...

Perform `exs` only on rank 0, otherwise know as "root" rank.
Other ranks will wait for the root rank to finish before continuing
"""
macro root(exp)
    command = quote
        if ClimaOcean.mpi_initialized()
            rank = ClimaOcean.mpi_rank()
            if rank == 0
                $exp
            end
            ClimaOcean.global_barrier()
        else
            $exp
        end
    end
    return esc(command)
end

"""
    @onrank rank, exs...

Perform `exp` only on rank `rank`
Other ranks will wait for the root rank to finish before continuing.
The expression is run anyways if MPI in not initialized
"""
macro onrank(exp_with_rank)
    on_rank = exp_with_rank.args[1]
    exp  = exp_with_rank.args[2]
    command = quote
        mpi_initialized = ClimaOcean.mpi_initialized()
        rank = ClimaOcean.mpi_rank()
        if !mpi_initialized
            $exp
        else
            if rank == $on_rank
                $exp
            end
            ClimaOcean.global_barrier()
        end
    end

    return esc(command)
end

""" 
    @distribute for i in iterable
        ...
    end

Distribute a `for` loop among different ranks
"""
macro distribute(exp)
    if exp.head != :for
        error("The `@distribute` macro expects a `for` loop")
    end

    iterable = exp.args[1].args[2]
    variable = exp.args[1].args[1]
    forbody  = exp.args[2]

    new_loop = quote
        mpi_initialized = ClimaOcean.mpi_initialized()
        if !mpi_initialized
            $exp
        else
            rank   = ClimaOcean.mpi_rank()
            nprocs = ClimaOcean.mpi_size()
            for (counter, $variable) in enumerate($iterable)
                if (counter - 1) % nprocs == rank
                    $forbody
                end
            end
            ClimaOcean.global_barrier()
        end
    end

    return esc(new_loop)
end

"""
    @handshake exs...

perform `exs` on all ranks, but only one rank at a time, where
ranks `r2 > r1` wait for rank `r1` to finish before executing `exs`
"""
macro handshake(exp)
    command = quote
        mpi_initialized = ClimaOcean.mpi_initialized()
        if !mpi_initialized
            $exp
        else
            rank   = ClimaOcean.mpi_rank()
            nprocs = ClimaOcean.mpi_size()
            for r in 0 : nprocs -1
                if rank == r
                    $exp
                end
                ClimaOcean.global_barrier()
            end
        end
    end
    return esc(command)
end