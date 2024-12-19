include("runtests_setup.jl")

using MPI
MPI.Init()

using NCDatasets
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

@testset "Distributed Bathymetry interpolation" begin
    # We start by building a fake bathyemtry on rank 0 and save it to file
    @root begin
        λ = 0:1:360
        φ = 0:1:20

        ds = NCDataset("./trivial_bathymetry.nc","c")

        # Define the dimension "lon" and "lat" with the size 361 and 21 resp.
        defDim(ds, "lon", 361)
        defDim(ds, "lat", 21)


        # Define the variables z
        z = defVar(ds, "z", Float32, ("lon","lat"))

        # Generate some example data
        data = [Float32(-i) for i = 1:361, j = 1:21]

        # write a the complete data set
        v[:,:] = data

        close(ds)
    end

    arch = Distributed(CPU(), partition=Partition(4, 1))

    grid = LatitudeLongitudeGrid(arch;
                                 size = (80, 10, 1),
                                 longitude = (0, 360),
                                 latitude = (0, 20),
                                 z = (0, 1))

    # It is linear interpolation so the data should be exactly equal to `i`
    bottom_height = regrid_bathymetry(grid; interpolation_passes=10)
    Nx, Ny, _ = size(grid)
    rank      = arch.local_rank
    
    irange = rank * Nx + 1 : (rank + 1) * Nx
    data   = repeat(irange, 1, Ny)
    @test interior(bottom_height) == - data
end