include("runtests_setup.jl")

@testset "JRA55 and data wrangling utilities" begin
    for arch in test_architectures

        A = typeof(arch)
        @info "Testing jra55_field_time_series on $A..."

        # This should download a file called "RYF.rsds.1990_1991.nc"
        time_indices = 1:3
        source_fts = ClimaOcean.JRA55.jra55_field_time_series(:shortwave_radiation, arch; time_indices)

        shortwave_radiation_filename = "RYF.rsds.1990_1991.nc"
        @test isfile(shortwave_radiation_filename)
        rm(shortwave_radiation_filename)

        @test source_fts isa FieldTimeSeries
        @test source_fts.grid isa LatitudeLongitudeGrid

        Nx, Ny, Nz, Nt = size(source_fts)
        @test Nx == 640
        @test Ny == 320
        @test Nz == 1
        @test Nt == length(time_indices)

        @info "Testing interpolate_field_time_series! on $A..."
        # Make target grid and field
        resolution = 1 # degree, eg 1/4
        Nxt = Int(360 / resolution)

        southern_limit = -79
        northern_limit = -30
        j₁ = (90 + southern_limit) / resolution
        j₂ = (90 + northern_limit) / resolution + 1
        Nyt = Int(j₂ - j₁ + 1)

        target_grid = LatitudeLongitudeGrid(arch,
                                            size = (Nxt, Nyt, 1);
                                            longitude = (0, 360),
                                            latitude = (southern_limit, northern_limit),
                                            z = (0, 1),
                                            topology = (Periodic, Bounded, Bounded))

        times = source_fts.times
        target_fts = FieldTimeSeries{Center, Center, Nothing}(target_grid, times)

        ClimaOcean.DataWrangling.interpolate_field_time_series!(target_fts, source_fts)
        # Random regression test
        CUDA.@allowscalar begin
            @test target_fts[1, 1, 1, 1] == 222.24310434509874    
        end

        @test target_fts.times == source_fts.times

        # What else might we test?

        @info "Testing save_field_time_series! on $A..."
        filepath = "JRA55_shortwave_radiation_test.jld2"
        ClimaOcean.DataWrangling.save_field_time_series!(target_fts, path=filepath, name="Qsw")
        @test isfile(filepath)

        # Test that we can load the data back
        Qswt = FieldTimeSeries(filepath, "Qsw")
        @test parent(Qswt) == parent(target_fts)    
        @test Qswt.times == target_fts.times
    end 
end

