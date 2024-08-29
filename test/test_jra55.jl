include("runtests_setup.jl")

using ClimaOcean.JRA55: download_jra55_cache
using ClimaOcean.OceanSeaIceModels: PrescribedAtmosphere

@testset "JRA55 and data wrangling utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing reanalysis_field_time_series on $A..."

        test_name = :downwelling_shortwave_radiation
        time_indices = 1:3

        # This should download a file called "RYF.rsds.1990_1991.nc"
        jra55_fts = JRA55_field_time_series(test_name; architecture=arch, time_indices)

        test_filename = joinpath(download_jra55_cache, "RYF.rsds.1990_1991.nc")

        @test jra55_fts isa FieldTimeSeries
        @test jra55_fts.grid isa LatitudeLongitudeGrid
        
        Nx, Ny, Nz, Nt = size(jra55_fts)
        @test Nx == 640
        @test Ny == 320
        @test Nz == 1
        @test Nt == length(time_indices)

        CUDA.@allowscalar begin
            @test jra55_fts[1, 1, 1, 1]   == 430.98105f0
            @test jra55_fts[641, 1, 1, 1] == 430.98105f0
        end

        # Test that halo regions were filled to respect boundary conditions
        CUDA.@allowscalar begin
            @test view(jra55_fts.data, 1, :, 1, :) == view(jra55_fts.data, Nx+1, :, 1, :)
        end

        @info "Testing loading preprocessed JRA55 data on $A..."
        in_memory_jra55_fts = JRA55_field_time_series(test_name;
                                                      time_indices,
                                                      architecture = arch,
                                                      backend = InMemory(2))

        @test in_memory_jra55_fts isa FieldTimeSeries

        @test interior(in_memory_jra55_fts[1]) == interior(jra55_fts[1])

        # Clean up
        rm(in_memory_jra55_fts.path)

        @info "Testing interpolate_field_time_series! on $A..."
        # Make target grid and field
        resolution = 1 # degree, eg 1/4
        Nx = Int(360 / resolution)

        southern_limit = -79
        northern_limit = -30
        j₁ = (90 + southern_limit) / resolution
        j₂ = (90 + northern_limit) / resolution + 1
        Ny = Int(j₂ - j₁ + 1)

        target_grid = LatitudeLongitudeGrid(arch,
                                            size = (Nx, Ny, 1);
                                            longitude = (0, 360),
                                            latitude = (southern_limit, northern_limit),
                                            z = (0, 1),
                                            topology = (Periodic, Bounded, Bounded))

        times = jra55_fts.times
        boundary_conditions = jra55_fts.boundary_conditions
        target_fts = FieldTimeSeries{Center, Center, Nothing}(target_grid, times; boundary_conditions)
        interpolate!(target_fts, jra55_fts)

        # Random regression test
        CUDA.@allowscalar begin
            @test target_fts[1, 1, 1, 1] == 222.243136478611

            # Only include this if we are filling halo regions within
            # interpolate_field_time_series
            @test target_fts[Nx + 1, 1, 1, 1] == 222.243136478611
        end

        @test target_fts.times == jra55_fts.times

        # What else might we test?

        @info "Testing save_field_time_series! on $A..."
        filepath = "JRA55_downwelling_shortwave_radiation_test.jld2"
        ClimaOcean.DataWrangling.save_field_time_series!(target_fts, path=filepath, name="Qsw",
                                                         overwrite_existing = true)
        @test isfile(filepath)

        # Test that we can load the data back
        Qswt = FieldTimeSeries(filepath, "Qsw")
        @test parent(Qswt.data) == parent(target_fts.data)
        @test Qswt.times == target_fts.times
        rm(filepath)

        #####
        ##### JRA55 prescribed atmosphere
        #####

        backend    = JRA55NetCDFBackend(2) 
        atmosphere = JRA55_prescribed_atmosphere(arch; backend, include_rivers_and_icebergs=false)
        @test atmosphere isa PrescribedAtmosphere
        @test isnothing(atmosphere.auxiliary_freshwater_flux)

        atmosphere = JRA55_prescribed_atmosphere(arch; backend, include_rivers_and_icebergs=true)
        @test haskey(atmosphere.auxiliary_freshwater_flux, :rivers)
        @test haskey(atmosphere.auxiliary_freshwater_flux, :icebergs)
    end 
end

