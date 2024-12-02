include("runtests_setup.jl")

using MPI
MPI.Init()

using ClimaOcean.ECCO: download_dataset, metadata_path
using CFTime
using Dates

@testset begin
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    @onrank 0, begin
        @test rank == 0
    end

    @root begin
        @test rank == 0
    end

    @onrank 1, begin
        @test rank == 1
    end

    @onrank 2, begin
        @test rank == 2
    end

    @onrank 3, begin
        @test rank == 3
    end

    a = Int[]
    
    @distribute for i in 1:10
        push!(a, i)
    end

    @root begin
        @test a == [1, 5, 9]
    end

    @onrank 1, begin
        @test a == [2, 6, 10]
    end

    @onrank 2, begin
        @test a == [3, 7]
    end

    @onrank 3, begin
        @test a == [4, 8]
    end
end

@testset "Distributed ECCO download" begin
    dates = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(1994, 4, 1)
    metadata = ECCOMetadata(:u_velocity; dates)
    download_dataset(metadata)

    @root for metadatum in metadata
        @test isfile(metadata_path(metadatum))
    end
end
