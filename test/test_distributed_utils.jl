include("runtests_setup.jl")

using MPI
MPI.Init()

using ClimaOcean.DataWrangling: metadata_path
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: reconstruct_global_grid
using CFTime
using Dates
using NCDatasets

# We start by building a fake bathymetry on rank 0 and save it to file
rm("./trivial_bathymetry.nc", force=true)

res = 0.5 # degrees
λ = -180+res/2:res:180-res/2
φ = 0:res:50

Nλ = length(λ)
Nφ = length(φ)

@root begin
    ds = NCDataset("./trivial_bathymetry.nc", "c")

    # Define the dimension "lon" and "lat" with the size Nλ and Nφ respectively
    defDim(ds, "lon", Nλ)
    defDim(ds, "lat", Nφ)
    defVar(ds, "lat", Float32, ("lat", ))
    defVar(ds, "lon", Float32, ("lon", ))

    # Define the variables z
    z = defVar(ds, "z", Float32, ("lon", "lat"))

    # Generate some example data
    data = [Float32(-i) for i = 1:Nλ, j = 1:Nφ]

    # write a the complete data set
    ds["lon"][:] = λ
    ds["lat"][:] = φ
    z[:, :] = data

    close(ds)
end

struct TrivalBathymetry end

import ClimaOcean.DataWrangling: download_dataset, z_interfaces, longitude_interfaces, latitude_interfaces, metadata_filename

download_dataset(::Metadatum{<:TrivalBathymetry, Nothing, Nothing}) = nothing
Base.size(::TrivalBathymetry) = (Nλ, Nφ, 1)
Base.size(::TrivalBathymetry, variable) = (Nλ, Nφ, 1)
z_interfaces(::TrivalBathymetry) = (0, 1)
longitude_interfaces(::TrivalBathymetry) = (-180, 180)
latitude_interfaces(::TrivalBathymetry) = (0, 50)
metadata_filename(metadatum::Metadatum{<:TrivalBathymetry, Nothing, Nothing}) = "trivial_bathymetry.nc"

@testset "Distributed ECCO download" begin
    dates = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(1994, 4, 1)
    metadata = Metadata(:u_velocity; dataset=ECCO4Monthly(), dates)
    download_dataset(metadata)

    @root for metadatum in metadata
        @test isfile(metadata_path(metadatum))
    end
end

@testset "Distributed Bathymetry interpolation" begin
    TrivialBathymetry_metadata = Metadata(:z, TrivalBathymetry(), nothing, nothing, ".")

    global_grid = LatitudeLongitudeGrid(CPU();
                                        size = (40, 40, 1),
                                        longitude = (0, 100),
                                        latitude = (0, 20),
                                        z = (0, 1))

    interpolation_passes = 4
    global_height = regrid_bathymetry(global_grid, TrivialBathymetry_metadata;
                                      interpolation_passes)

    arch_x  = Distributed(CPU(), partition=Partition(4, 1))
    arch_y  = Distributed(CPU(), partition=Partition(1, 4))
    arch_xy = Distributed(CPU(), partition=Partition(2, 2))

    for arch in (arch_x, arch_y, arch_xy)
        local_grid = LatitudeLongitudeGrid(arch;
                                           size = (40, 40, 1),
                                           longitude = (0, 100),
                                           latitude = (0, 20),
                                           z = (0, 1))

        local_height = regrid_bathymetry(local_grid, TrivialBathymetry_metadata;
                                         interpolation_passes)

        Nx, Ny, _ = size(local_grid)
        rx, ry, _ = arch.local_index
        irange    = (rx - 1) * Nx + 1 : rx * Nx
        jrange    = (ry - 1) * Ny + 1 : ry * Ny

        begin
            @test interior(global_height, irange, jrange, 1) == interior(local_height, :, :, 1)
        end
    end
end

@testset "Distributed atmospheric halo filling" begin
    function create_sample_output(arch, filename)

        Nx, Ny, Nz = 180, 100, 20

        z = ExponentialDiscretization(Nz, -6000, 0, mutable=true)
        underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), z, halo = (7, 7, 4))

        ETOPOmetadata = Metadatum(:bottom_height, dataset=ETOPO2022())
        ClimaOcean.DataWrangling.download_dataset(ETOPOmetadata)

        bottom_height = regrid_bathymetry(underlying_grid, ETOPOmetadata;
                                          minimum_depth = 15,
                                          interpolation_passes = 1,
                                          major_basins = 1)

        grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map=true)

        ocean = ocean_simulation(grid; free_surface = SplitExplicitFreeSurface(grid; substeps = 70))

        radiation  = Radiation(arch)
        atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(100), include_rivers_and_icebergs=true)

        coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

        simulation = Simulation(coupled_model; Δt=60, stop_time=11minutes)

        dates = vcat(collect(DateTime(1991, 1, 1): Month(1): DateTime(1991, 5, 1)),
                     collect(DateTime(1990, 5, 1): Month(1): DateTime(1990, 12, 1)))

        dataset = EN4Monthly()

        temperature = Metadata(:temperature; dates, dataset = dataset)
        salinity    = Metadata(:salinity;    dates, dataset = dataset)

        set!(ocean.model, T=Metadata(:temperature; dates=first(dates), dataset = dataset),
                          S=Metadata(:salinity;    dates=first(dates), dataset = dataset))

        ## Print a progress message
        progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                        iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                        maximum(abs, sim.model.ocean.model.velocities.w), prettytime(sim.run_wall_time))

        add_callback!(simulation, progress_message, IterationInterval(1))

        tracers = ocean.model.tracers
        velocities = ocean.model.velocities

        outputs = merge(tracers, velocities)

        simulation.output_writers[:snapshot] = JLD2Writer(ocean.model, outputs;
                                                          schedule = TimeInterval(10minutes),
                                                          filename,
                                                          indices = (:, :, Nz),
                                                          with_halos = false,
                                                          overwrite_existing = true,
                                                          array_type = Array{Float32})

        run!(simulation)
    end

    
    filename1 = "snapshot_distributed"
    arch = Distributed(CPU(); partition = Partition(y = DistributedComputations.Equal()), synchronized_communication=true)
    create_sample_output(arch, filename1)

    filename2 = "snapshot_serial"
    arch = CPU()
    create_sample_output(arch, filename2)

    distributed_files = filter(f -> occursin("_rank", f), glob("$filename1*.jld2"))

    serial_file = glob("$filename2.jld2")

    ranks = size(distributed_files)
    var = "T"
    T_rank_dist = []

    for rank in 0:ranks-1
        fname_rank = "snapshot_distributed_rank$(rank).jld2"
        @info "Reconstructing global grid from $fname_rank"
        keys_iters= keys(jldopen(fname_rank)["timeseries"][var])[2:end]
        T_rank_full = jldopen(fname_rank)["timeseries"][var][keys_iters[lastindex(keys_iters)]][:,:,1]
        push!(T_rank_dist, T_rank_full)
    end
    T_rank_dist_all = cat(T_rank_dist...; dims = 2)

    field2 = jldopen(serial_file[1])
    timesteps = keys(jldopen(serial_file[1])["timeseries"][var])[end]
    T_serial = field2["timeseries"][var][timesteps][:,:,1]

    @test (maximum(abs.(T_serial .- T_rank_dist_all)) < 1e-10)
end

MPI.Finalize()
