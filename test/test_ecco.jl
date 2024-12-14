include("runtests_setup.jl")

using CFTime
using Dates
using ClimaOcean

using ClimaOcean.ECCO
using ClimaOcean.ECCO: ECCO_field, metadata_path, ECCO_times
using ClimaOcean.DataWrangling: NearestNeighborInpainting

using Oceananigans.Grids: topology
using Oceananigans.OutputReaders: time_indices
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.Units

using CUDA: @allowscalar

start_date = DateTimeProlepticGregorian(1993, 1, 1)
end_date = DateTimeProlepticGregorian(1993, 2, 1)
dates = start_date : Month(1) : end_date

# Inpaint only the first two cells inside the missing mask
inpainting = NearestNeighborInpainting(2)

@testset "ECCO fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        for name in (:temperature, :salinity)
            @info "Testing ECCO_field on $A..."
            metadata = ECCOMetadata(name, dates, ECCO4Monthly())
            restoring = ECCORestoring(metadata; rate = 1 / 1000.0, inpainting)

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

            @test fts.times[1] == ECCO_times(metadata)[1]
            @test fts.times[end] == ECCO_times(metadata)[end]

            datum = first(metadata)
            ψ = ECCO_field(datum, architecture=arch, inpainting=NearestNeighborInpainting(2))
            datapath = ClimaOcean.DataWrangling.ECCO.inpainted_metadata_path(datum)
            @test isfile(datapath)
        end
    end
end

@testset "Inpainting algorithm" begin
    for arch in test_architectures
        T_metadata = ECCOMetadata(:temperature, dates[1], ECCO4Monthly())

        grid = LatitudeLongitudeGrid(arch,
                                     size = (100, 100, 10),
                                     latitude = (-75, 75),
                                     longitude = (0, 360),
                                     z = (-200, 0),
                                     halo = (6, 6, 6))

        fully_inpainted_field = CenterField(grid)
        partially_inpainted_field = CenterField(grid)

        set!(fully_inpainted_field,     T_metadata; inpainting = NearestNeighborInpainting(Inf))
        set!(partially_inpainted_field, T_metadata; inpainting = NearestNeighborInpainting(1))

        fully_inpainted_interior = on_architecture(CPU(), interior(fully_inpainted_field))
        partially_inpainted_interior = on_architecture(CPU(), interior(partially_inpainted_field))

        @test all(fully_inpainted_interior .!= 0)
        @test any(partially_inpainted_interior .== 0)
    end
end

@testset "LinearlyTaperedPolarMask" begin
    for arch in test_architectures
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

        t_restoring = ECCORestoring(:temperature, arch;
                                    dates,
                                    mask,
                                    rate = 1 / 1000.0,
                                    inpainting)

        fill!(t_restoring.field_time_series[1], 1.0)
        fill!(t_restoring.field_time_series[2], 1.0)

        T = CenterField(grid)
        fields = (; T)
        clock  = Clock(; time = 0)

        @test @allowscalar t_restoring(1, 1,   10, grid, clock, fields) == t_restoring.rate
        @test @allowscalar t_restoring(1, 11,  10, grid, clock, fields) == t_restoring.rate / 2
        @test @allowscalar t_restoring(1, 21,  10, grid, clock, fields) == 0
        @test @allowscalar t_restoring(1, 80,  10, grid, clock, fields) == 0
        @test @allowscalar t_restoring(1, 90,  10, grid, clock, fields) == t_restoring.rate / 2
        @test @allowscalar t_restoring(1, 100, 10, grid, clock, fields) == t_restoring.rate
        @test @allowscalar t_restoring(1, 1,   5,  grid, clock, fields) == 0
        @test @allowscalar t_restoring(1, 10,  5,  grid, clock, fields) == 0
    end
end

@testset "Setting a field with ECCO" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size=(10, 10, 10),
                                     latitude=(-60, -40),
                                     longitude=(10, 15), z=(-200, 0))

        field = CenterField(grid)

        @test begin
            set!(field, ECCOMetadata(:temperature))
            set!(field, ECCOMetadata(:salinity))
            true
        end
    end
end

@testset "Timestepping with ECCORestoring" begin
    for arch in test_architectures

        grid  = LatitudeLongitudeGrid(arch;
                                      size = (10, 10, 10),
                                      latitude = (-60, -40),
                                      longitude = (10, 15),
                                      z = (-200, 0),
                                      halo = (6, 6, 6))

        field = CenterField(grid)

        @test begin
            set!(field, ECCOMetadata(:temperature))
            set!(field, ECCOMetadata(:salinity))
            true
        end

        forcing_T = ECCORestoring(:temperature, arch;
                                  dates,
                                  rate = 1 / 1000.0,
                                  inpainting)

        ocean = ocean_simulation(grid; forcing = (; T = forcing_T))

        @test begin
            time_step!(ocean)
            time_step!(ocean)
            true
        end
    end
end

@testset "Setting temperature and salinity to ECCO" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch; 
                                     size=(10, 10, 10),
                                     latitude=(-60, -40),
                                     longitude=(10, 15),
                                     z=(-200, 0),
                                     halo = (7, 7, 7))

        ocean = ocean_simulation(grid)
        date = DateTimeProlepticGregorian(1993, 1, 1)
        set!(ocean.model, T=ECCOMetadata(:temperature, date), S=ECCOMetadata(:salinity, date))
    end
end

@testset "ECCO dataset cycling boundaries" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size=(10, 10, 10),
                                     latitude=(-60, -40),
                                     longitude=(10, 15),
                                     z=(-200, 0),
                                     halo = (7, 7, 7))

        start_date = DateTimeProlepticGregorian(1993, 1, 1)
        end_date = DateTimeProlepticGregorian(1993, 5, 1)
        dates = start_date : Month(1) : end_date

        t_restoring = ECCORestoring(arch, :temperature;
                                    dates,
                                    rate = 1 / 1000.0,
                                    inpainting)

        times = ECCO_times(t_restoring.field_time_series.backend.metadata)
        ocean = ocean_simulation(grid, forcing = (; T = t_restoring))

        ocean.model.clock.time = times[3] + 2 * Units.days
        update_state!(ocean.model)

        @test t_restoring.field_time_series.backend.start == 3

        # Compile
        time_step!(ocean)

        # Try stepping out of the ECCO dataset bounds
        ocean.model.clock.time = last(times) + 2 * Units.days

        update_state!(ocean.model)
        
        @test begin
            time_step!(ocean)
            true
        end

        # The backend has cycled to the end
        @test time_indices(t_restoring.field_time_series) == (5, 1)
    end
end
