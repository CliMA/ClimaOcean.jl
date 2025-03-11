include("runtests_setup.jl")

using MPI
MPI.Init()

using NCDatasets
using ClimaOcean.ECCO: download_dataset, metadata_path
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: reconstruct_global_grid
using CFTime
using Dates

@testset "Distributed ECCO download" begin
    dates = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(1994, 4, 1)
    metadata = Metadata(:u_velocity; dataset=ECCO4Monthly(), dates)
    download_dataset(metadata)

    @root for metadatum in metadata
        @test isfile(metadata_path(metadatum))
    end
end

@testset "Distributed Bathymetry interpolation" begin
    # We start by building a fake bathymetry on rank 0 and save it to file
    @root begin
        λ = -179.5:0.1:179.5
        φ = 0:0.1:50

        Nλ = length(λ)
        Nφ = length(φ)

        ds = NCDataset("./trivial_bathymetry.nc", "c")

        # Define the dimension "lon" and "lat" with the size 361 and 51 resp.
        defDim(ds, "lon", Nλ)
        defDim(ds, "lat", Nφ)
        defVar(ds, "lat", Float32, ("lat", ))
        defVar(ds, "lon", Float32, ("lon", ))

        # Define the variables z
        z = defVar(ds, "z", Float32, ("lon","lat"))

        # Generate some example data
        data = [Float32(-i) for i = 1:Nλ, j = 1:Nφ]

        # write a the complete data set
        ds["lon"][:] = λ
        ds["lat"][:] = φ
        z[:,:] = data
        
        close(ds)
    end

    global_grid = LatitudeLongitudeGrid(CPU();
                                        size = (40, 40, 1),
                                        longitude = (0, 100),
                                        latitude = (0, 20),
                                        z = (0, 1))

    global_height = regrid_bathymetry(global_grid; 
                                    dir = "./",
                                    filename = "trivial_bathymetry.nc",
                                    interpolation_passes=10)

    arch_x  = Distributed(CPU(), partition=Partition(4, 1))
    arch_y  = Distributed(CPU(), partition=Partition(1, 4))
    arch_xy = Distributed(CPU(), partition=Partition(2, 2))

    for arch in (arch_x, arch_y, arch_xy)
        local_grid = LatitudeLongitudeGrid(arch;
                                           size = (40, 40, 1),
                                           longitude = (0, 100),
                                           latitude = (0, 20),
                                            z = (0, 1))

        local_height = regrid_bathymetry(local_grid; 
                                         dir = "./",
                                         filename = "trivial_bathymetry.nc",
                                         interpolation_passes=10)

        Nx, Ny, _ = size(local_grid)
        rx, ry, _ = arch.local_index
        irange    = (rx - 1) * Nx + 1 : rx * Nx
        jrange    = (ry - 1) * Ny + 1 : ry * Ny
        
        @handshake begin
            @test interior(global_height, irange, jrange, 1) == interior(local_height, :, :, 1)
        end
    end
end