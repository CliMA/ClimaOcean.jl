include("runtests_setup.jl")

using MPI
MPI.Init()

using ClimaOcean.ECCO: download_dataset, metadata_path
using CFTime
using Dates

@testset begin
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    @onrank 0 begin
        @test rank == 0
    end

    @root begin
        @test rank == 0
    end

    @onrank 1 begin
        @test rank == 1
    end

    @onrank 2 begin
        @test rank == 2
    end

    @onrank 3 begin
        @test rank == 3
    end

    a = Int[]
    
    @distribute for i in 1:10
        push!(a, i)
    end

    @root begin
        @test a == [1, 5, 9]
    end

    @onrank 1 begin
        @test a == [2, 6, 10]
    end

    @onrank 2 begin
        @test a == [3, 7]
    end

    @onrank 3 begin
        @test a == [4, 8]
    end
  

    split_comm = MPI.Comm_split(MPI.COMM_WORLD, rank % 2, rank)

    a = Int[]
    
    @distribute split_comm for i in 1:10
        push!(a, i)
    end

    @onrank split_comm 0 @test a == [1, 3, 5, 7, 9]
    @onrank split_comm 1 @test a == [2, 4, 6, 8, 10]
end

@testset "Distributed ECCO download" begin
    dates = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(1994, 4, 1)
    metadata = ECCOMetadata(:u_velocity; dates)
    download_dataset(metadata)

    @root for metadatum in metadata
        @test isfile(metadata_path(metadatum))
    end
end
