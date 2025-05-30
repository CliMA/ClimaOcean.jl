using ClimaOcean
using Oceananigans
using CUDA
using Test

using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: metadata_path
using ClimaOcean.EN4
using ClimaOcean.ECCO
using ClimaOcean.ETOPO
using ClimaOcean.JRA55

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.OutputReaders: interpolate!

using CFTime
using Dates

using CUDA: @allowscalar

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]

start_date = DateTimeProlepticGregorian(1993, 1, 1)

test_datasets = (ECCO2Monthly(), ECCO2Daily(), ECCO4Monthly(), EN4Monthly())

#####
##### Test utilities
#####

function test_setting_from_metadata(arch, dataset, start_date, inpainting,
                                    loc = (Center, Center, Center),
                                    names = (:temperature, :salinity))

    grid = LatitudeLongitudeGrid(arch;
                                 size = (10, 10, 10),
                                 latitude = (-60, -40),
                                 longitude = (10, 15),
                                 z = (-200, 0))

    field = Field(loc, grid)

    @test begin
        for name in names
            set!(field, Metadatum(name; dataset, date=start_date); inpainting)
        end
        true
    end

    return nothing
end

function test_timestepping_with_dataset(arch, dataset, start_date, inpainting)
    grid  = LatitudeLongitudeGrid(arch;
                                  size = (10, 10, 10),
                                  latitude = (-60, -40),
                                  longitude = (10, 15),
                                  z = (-200, 0),
                                  halo = (6, 6, 6))

    field = CenterField(grid)

    @test begin
        set!(field, Metadatum(:temperature; dataset, date=start_date); inpainting)
        set!(field, Metadatum(:salinity; dataset, date=start_date); inpainting)
        true
    end

    ocean = ocean_simulation(grid; verbose=false)

    @test begin
        time_step!(ocean)
        time_step!(ocean)
        true
    end

    return nothing
end

function test_ocean_metadata_utilities(arch, dataset, dates, inpainting)
    for name in (:temperature, :salinity)
        metadata = Metadata(name; dates, dataset)
        restoring = DatasetRestoring(metadata, arch; rate=1/1000, inpainting)

        for datum in metadata
            @test isfile(metadata_path(datum))
        end

        fts = restoring.field_time_series
        @test fts isa FieldTimeSeries
        @test fts.grid isa LatitudeLongitudeGrid
        @test topology(fts.grid) == (Periodic, Bounded, Bounded)

        Nx, Ny, Nz = size(interior(fts))
        Nt = length(fts.times)

        @test Nx == size(metadata)[1]
        @test Ny == size(metadata)[2]
        @test Nz == size(metadata)[3]
        @test Nt == size(metadata)[4]

        @test @allowscalar fts.times[1] == native_times(metadata)[1]
        @test @allowscalar fts.times[end] == native_times(metadata)[end]

        datum = first(metadata)
        ψ = Field(datum, arch, inpainting=NearestNeighborInpainting(2))
        datapath = ClimaOcean.DataWrangling.inpainted_metadata_path(datum)
        @test isfile(datapath)
    end

    return nothing
end

function test_dataset_restoring(arch, dataset, dates, inpainting)
    grid = LatitudeLongitudeGrid(arch;
                                 size = (100, 100, 10),
                                 latitude = (-75, 75),
                                 longitude = (0, 360),
                                 z = (-200, 0),
                                 halo = (6, 6, 6))

    φ₁ = @allowscalar grid.φᵃᶜᵃ[1]
    φ₂ = @allowscalar grid.φᵃᶜᵃ[21]
    φ₃ = @allowscalar grid.φᵃᶜᵃ[80]
    φ₄ = @allowscalar grid.φᵃᶜᵃ[100]
    z₁ = @allowscalar grid.z.cᵃᵃᶜ[6]

    mask = LinearlyTaperedPolarMask(northern = (φ₃, φ₄),
                                    southern = (φ₁, φ₂),
                                    z = (z₁, 0))

    for name in (:temperature, :salinity)
        metadata = Metadata(name; dates, dataset)
        var_restoring = DatasetRestoring(metadata, arch; mask, inpainting, rate=1/1000)

        fill!(var_restoring.field_time_series[1], 1.0)
        fill!(var_restoring.field_time_series[2], 1.0)

        T = CenterField(grid)
        S = CenterField(grid)
        fields = (; T, S)
        clock  = Clock(; time = 0)

        @allowscalar begin
            @test var_restoring(1, 1,   10, grid, clock, fields) == var_restoring.rate
            @test var_restoring(1, 11,  10, grid, clock, fields) == var_restoring.rate / 2
            @test var_restoring(1, 21,  10, grid, clock, fields) == 0
            @test var_restoring(1, 80,  10, grid, clock, fields) == 0
            @test var_restoring(1, 90,  10, grid, clock, fields) == var_restoring.rate / 2
            @test var_restoring(1, 100, 10, grid, clock, fields) == var_restoring.rate
            @test var_restoring(1, 1,   5,  grid, clock, fields) == 0
            @test var_restoring(1, 10,  5,  grid, clock, fields) == 0
        end
    end

    return nothing
end

function test_timestepping_with_dataset_restoring(arch, dataset, dates, inpainting)
    grid = LatitudeLongitudeGrid(arch;
                                 size = (10, 10, 10),
                                 latitude = (-60, -40),
                                 longitude = (10, 15),
                                 z = (-200, 0),
                                 halo = (6, 6, 6))

    T_meta = Metadata(:temperature; dates, dataset)
    T_forcing = DatasetRestoring(T_meta, arch; inpainting, rate=1/1000)

    ocean = ocean_simulation(grid; forcing = (; T = T_forcing), verbose=false)

    @test begin
        time_step!(ocean)
        time_step!(ocean)
        true
    end

    return nothing
end

function test_cycling_dataset_restoring(arch, dataset, dates, inpainting)
    grid = LatitudeLongitudeGrid(arch;
                                 size = (10, 10, 10),
                                 latitude = (-60, -40),
                                 longitude = (10, 15),
                                 z = (-200, 0),
                                 halo = (7, 7, 7))

    time_indices_in_memory = 2
    start_date = dates[1]
    end_date = dates[end]

    Tmetadata1 = Metadata(:temperature; dates, dataset)
    Tmetadata2 = Metadata(:temperature; start_date, end_date, dataset)

    for Tmetadata in (Tmetadata1, Tmetadata2)
        T_restoring = DatasetRestoring(Tmetadata, arch; time_indices_in_memory, inpainting, rate=1/1000)

        times = native_times(T_restoring.field_time_series.backend.metadata)
        ocean = ocean_simulation(grid, forcing = (; T = T_restoring))

        # start a bit after time_index
        time_index = 3
        time_interval = dataset isa ECCO2Daily ? Units.hour : Units.day
        ocean.model.clock.time = times[time_index] + 2 * time_interval
        update_state!(ocean.model)

        @test time_indices(T_restoring.field_time_series) ==
            Tuple(range(time_index, length=time_indices_in_memory))

        @test T_restoring.field_time_series.backend.start == time_index

        # Compile
        time_step!(ocean)

        # Try stepping out of the dataset bounds
        # start a bit after last time_index
        ocean.model.clock.time = last(times) + 2 * time_interval

        update_state!(ocean.model)

        @test begin
            time_step!(ocean)
            true
        end

        # The backend has cycled to the end
        @test time_indices(T_restoring.field_time_series) ==
            mod1.(Tuple(range(length(times), length=time_indices_in_memory)), length(times))
    end
end

function test_inpainting_algorithm(arch, dataset, start_date, inpainting)
    T_metadatum = Metadatum(:temperature; dataset, date=start_date)

    grid = LatitudeLongitudeGrid(arch,
                                 size = (20, 20, 10),
                                 latitude = (-75, 75),
                                 longitude = (0, 360),
                                 z = (-4000, 0),
                                 halo = (6, 6, 6))

    fully_inpainted_field = CenterField(grid)
    partially_inpainted_field = CenterField(grid)

    set!(fully_inpainted_field, T_metadatum; inpainting = NearestNeighborInpainting(Inf))
    set!(partially_inpainted_field, T_metadatum; inpainting = NearestNeighborInpainting(1))

    fully_inpainted_interior = on_architecture(CPU(), interior(fully_inpainted_field))
    partially_inpainted_interior = on_architecture(CPU(), interior(partially_inpainted_field))

    @test all(fully_inpainted_interior .!= 0)
    @test any(partially_inpainted_interior .== 0)

    return nothing
end
