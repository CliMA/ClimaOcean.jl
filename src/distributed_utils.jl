using MPI

#####
##### Handle commands, typically downloading files
##### which should be executed by only one rank or distributed among ranks
#####

# Utilities to make the macro work importing only ClimaOcean and not MPI
mpi_initialized()    = MPI.Initialized()
mpi_rank(comm)       = MPI.Comm_rank(comm)
mpi_size(comm)       = MPI.Comm_size(comm)
global_barrier(comm) = MPI.Barrier(comm) 

"""
    @root communicator exs...

Perform `exs` only on rank 0 in communicator, otherwise known as the "root" rank.
Other ranks will wait for the root rank to finish before continuing.
If `communicator` is not provided, `MPI.COMM_WORLD` is used.
"""
macro root(communicator, exp)
    command = quote
        if ClimaOcean.mpi_initialized()
            rank = ClimaOcean.mpi_rank($communicator)
            if rank == 0
                $exp
            end
            ClimaOcean.global_barrier($communicator)
        else
            $exp
        end
    end
    return esc(command)
end

macro root(exp)
    command = quote
        @root MPI.COMM_WORLD $exp
    end
    return esc(command)
end

"""
    @onrank communicator rank exs...

Perform `exp` only on rank `rank` (0-based index) in `communicator`.
Other ranks will wait for rank `rank` to finish before continuing.
The expression is run anyways if MPI in not initialized.
If `communicator` is not provided, `MPI.COMM_WORLD` is used.
"""
macro onrank(communicator, on_rank, exp)
    command = quote
        mpi_initialized = ClimaOcean.mpi_initialized()
        if !mpi_initialized
            $exp
        else
            rank = ClimaOcean.mpi_rank($communicator)
            if rank == $on_rank
                $exp
            end
            ClimaOcean.global_barrier($communicator)
        end
    end

    return esc(command)
end

macro onrank(rank, exp)
    command = quote
        @onrank MPI.COMM_WORLD $rank $exp
    end
    return esc(command)
end

""" 
    @distribute communicator for i in iterable
        ...
    end

Distribute a `for` loop among different ranks in `communicator`.
If `communicator` is not provided, `MPI.COMM_WORLD` is used.
"""
macro distribute(communicator, exp)
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
            rank   = ClimaOcean.mpi_rank($communicator)
            nprocs = ClimaOcean.mpi_size($communicator)
            for (counter, $variable) in enumerate($iterable)
                if (counter - 1) % nprocs == rank
                    $forbody
                end
            end
            ClimaOcean.global_barrier($communicator)
        end
    end

    return esc(new_loop)
end

macro distribute(exp)
    command = quote
        @distribute MPI.COMM_WORLD $exp
    end
    return esc(command)
end

"""
    @handshake communicator exs...

perform `exs` on all ranks in `communicator`, but only one rank at a time, where
ranks `r2 > r1` wait for rank `r1` to finish before executing `exs`.
If `communicator` is not provided, `MPI.COMM_WORLD` is used.
"""
macro handshake(communicator, exp)
    command = quote
        mpi_initialized = ClimaOcean.mpi_initialized()
        if !mpi_initialized
            $exp
        else
            rank   = ClimaOcean.mpi_rank($communicator)
            nprocs = ClimaOcean.mpi_size($communicator)
            for r in 0 : nprocs -1
                if rank == r
                    $exp
                end
                ClimaOcean.global_barrier($communicator)
            end
        end
    end
    return esc(command)
end

macro handshake(exp)
    command = quote
        @handshake MPI.COMM_WORLD $exp
    end
    return esc(command)
end
