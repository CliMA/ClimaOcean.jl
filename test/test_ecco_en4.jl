for arch in test_architectures, dataset in test_datasets

    A = typeof(arch)
    @info "Testing $(typeof(dataset)) on $A"

    start_date = DateTime(1993, 1, 1)
    time_resolution = dataset isa ECCO2Daily ? Day(1) : Month(1)
    end_date = start_date + 4 * time_resolution
    dates = start_date : time_resolution : end_date

    @testset "Fields utilities" begin
        for name in (:temperature, :salinity)
            metadata = Metadata(name; dates, dataset)

            download_dataset(metadata) # just in case is not downloaded
            for datum in metadata
                @test isfile(metadata_path(datum))
            end

            datum = first(metadata)
            ψ = Field(datum, architecture=arch, inpainting=NearestNeighborInpainting(2))
            @test ψ isa Field
            datapath = ClimaOcean.DataWrangling.inpainted_metadata_path(datum)
            @test isfile(datapath)
        end
    end

    @testset "Inpainting algorithm" begin
        T_metadatum = Metadatum(:temperature; dataset, date=start_date)

        grid = LatitudeLongitudeGrid(arch,
                                     size = (100, 100, 10),
                                     latitude = (-75, 75),
                                     longitude = (0, 360),
                                     z = (-200, 0),
                                     halo = (6, 6, 6))

        fully_inpainted_field = CenterField(grid)
        partially_inpainted_field = CenterField(grid)

        set!(fully_inpainted_field,     T_metadatum; inpainting = NearestNeighborInpainting(Inf))
        set!(partially_inpainted_field, T_metadatum; inpainting = NearestNeighborInpainting(1))

        fully_inpainted_interior = on_architecture(CPU(), interior(fully_inpainted_field))
        partially_inpainted_interior = on_architecture(CPU(), interior(partially_inpainted_field))

        @test all(fully_inpainted_interior .!= 0)
        @test any(partially_inpainted_interior .== 0)
    end

    @testset "Setting a field from a dataset" begin
        grid = LatitudeLongitudeGrid(arch;
                                     size=(10, 10, 10),
                                     latitude=(-60, -40),
                                     longitude=(10, 15), z=(-200, 0))

        field = CenterField(grid)

        @test begin
            set!(field, Metadatum(:temperature; dataset, date=start_date))
            set!(field, Metadatum(:salinity;    dataset, date=start_date))
            true
        end
    end

    @testset "Setting temperature and salinity from dataset" begin
        grid = LatitudeLongitudeGrid(arch;
                                     size = (10, 10, 10),
                                     latitude = (-60, -40),
                                     longitude = (10, 15),
                                     z = (-200, 0),
                                     halo = (7, 7, 7))

        ocean = ocean_simulation(grid)

        set!(ocean.model, T=Metadatum(:temperature; dataset, date=start_date),
                          S=Metadatum(:salinity;    dataset, date=start_date))
    end

    @testset "Timestepping with fields from Dataset" begin
        grid  = LatitudeLongitudeGrid(arch;
                                      size = (10, 10, 10),
                                      latitude = (-60, -40),
                                      longitude = (10, 15),
                                      z = (-200, 0),
                                      halo = (6, 6, 6))

        field = CenterField(grid)

        @test begin
            set!(field, Metadatum(:temperature; dataset, date=start_date))
            set!(field, Metadatum(:salinity;    dataset, date=start_date))
            true
        end

        ocean = ocean_simulation(grid; verbose=false)

        @test begin
            time_step!(ocean)
            time_step!(ocean)
            true
        end
    end

    @testset "Field utilities" begin
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

            @test fts.times[1] == native_times(metadata)[1]
            @test fts.times[end] == native_times(metadata)[end]

            datum = first(metadata)
            ψ = Field(datum, architecture=arch, inpainting=NearestNeighborInpainting(2))
            datapath = ClimaOcean.DataWrangling.inpainted_metadata_path(datum)
            @test isfile(datapath)
        end
    end

    @testset "DatasetRestoring with LinearlyTaperedPolarMask" begin
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
            var_restoring = DatasetRestoring(name, arch; dataset, start_date, end_date, mask, inpainting, rate=1/1000)

            fill!(var_restoring.field_time_series[1], 1.0)
            fill!(var_restoring.field_time_series[2], 1.0)

            T = CenterField(grid)
            S = CenterField(grid)
            fields = (; T, S)
            clock  = Clock(; time = 0)

            @test @allowscalar var_restoring(1, 1,   10, grid, clock, fields) == var_restoring.rate
            @test @allowscalar var_restoring(1, 11,  10, grid, clock, fields) == var_restoring.rate / 2
            @test @allowscalar var_restoring(1, 21,  10, grid, clock, fields) == 0
            @test @allowscalar var_restoring(1, 80,  10, grid, clock, fields) == 0
            @test @allowscalar var_restoring(1, 90,  10, grid, clock, fields) == var_restoring.rate / 2
            @test @allowscalar var_restoring(1, 100, 10, grid, clock, fields) == var_restoring.rate
            @test @allowscalar var_restoring(1, 1,   5,  grid, clock, fields) == 0
            @test @allowscalar var_restoring(1, 10,  5,  grid, clock, fields) == 0
        end
    end

    @testset "Timestepping with DatasetRestoring" begin
        grid = LatitudeLongitudeGrid(arch;
                                     size = (10, 10, 10),
                                     latitude = (-60, -40),
                                     longitude = (10, 15),
                                     z = (-200, 0),
                                     halo = (6, 6, 6))

        field = CenterField(grid)

        @test begin
            set!(field, Metadatum(:temperature; dataset, date=start_date))
            set!(field, Metadatum(:salinity;    dataset, date=start_date))
            true
        end

        forcing_T = DatasetRestoring(:temperature, arch; dataset, start_date, end_date, inpainting, rate=1/1000)

        ocean = ocean_simulation(grid; forcing = (; T = forcing_T), verbose=false)

        @test begin
            time_step!(ocean)
            time_step!(ocean)
            true
        end
    end

    @testset "Dataset cycling boundaries" begin
        grid = LatitudeLongitudeGrid(arch;
                                     size = (10, 10, 10),
                                     latitude = (-60, -40),
                                     longitude = (10, 15),
                                     z = (-200, 0),
                                     halo = (7, 7, 7))

        start_date = DateTime(1993, 1, 1)
        time_resolution = dataset isa ECCO2Daily ? Day(1) : Month(1)
        end_date = start_date + 4 * time_resolution

        time_indices_in_memory = 2

        # create a T_restoring using both constructors
        T_restoring1 = DatasetRestoring(:temperature, arch; dataset, start_date, end_date, inpainting, rate=1/1000)

        Tmetadata = Metadata(:temperature; dates=start_date:time_resolution:end_date, dataset)
        T_restoring2 = DatasetRestoring(Tmetadata, arch; time_indices_in_memory, inpainting, rate=1/1000)

        for T_restoring in (T_restoring1, T_restoring2)
            times = native_times(T_restoring.field_time_series.backend.metadata)
            ocean = ocean_simulation(grid, forcing = (; T = T_restoring))

            time_interval = dataset isa ECCO2Daily ? 2 * Units.hours : 2 * Units.days

            ocean.model.clock.time = times[3] + time_interval
            update_state!(ocean.model)

            @test T_restoring.field_time_series.backend.start == 3

            # Compile
            time_step!(ocean)

            # Try stepping out of the dataset bounds
            ocean.model.clock.time = last(times) + time_interval

            update_state!(ocean.model)

            @test begin
                time_step!(ocean)
                true
            end

            # The backend has cycled to the end
            @test Oceananigans.OutputReaders.time_indices(T_restoring.field_time_series) == (5, 1)
        end
    end
end
