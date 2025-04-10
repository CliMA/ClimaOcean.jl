include("runtests_setup.jl")

using ClimaOcean.JRA55
using ClimaOcean.JRA55: download_JRA55_cache
using ClimaOcean.OceanSeaIceModels: PrescribedAtmosphere

@testset "JRA55 and data wrangling utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing reanalysis_field_time_series on $A..."

        # This should download files called "RYF.rsds.1990_1991.nc" and "RYF.tas.1990_1991.nc"
        for test_name in (:downwelling_shortwave_radiation, :temperature)
            dates = ClimaOcean.DataWrangling.all_dates(JRA55.RepeatYearJRA55(), test_name)
            end_date = dates[3]

            JRA55_fts = JRA55FieldTimeSeries(test_name, arch; end_date)
            test_filename = joinpath(download_JRA55_cache, "RYF.rsds.1990_1991.nc")

            @test JRA55_fts isa FieldTimeSeries
            @test JRA55_fts.grid isa LatitudeLongitudeGrid

            Nx, Ny, Nz, Nt = size(JRA55_fts)
            @test Nx == 640
            @test Ny == 320
            @test Nz == 1
            @test Nt == 3

            if test_name == :downwelling_shortwave_radiation
                CUDA.@allowscalar begin
                    @test JRA55_fts[1, 1, 1, 1]   == 430.98105f0
                    @test JRA55_fts[641, 1, 1, 1] == 430.98105f0
                end
            end

            # Test that halo regions were filled to respect boundary conditions
            CUDA.@allowscalar begin
                @test view(JRA55_fts.data, 1, :, 1, :) == view(JRA55_fts.data, Nx+1, :, 1, :)
            end

            @info "Testing Cyclical time_indices for JRA55 data on $A..."
            Nb = 4
            backend = JRA55NetCDFBackend(Nb)
            netcdf_JRA55_fts = JRA55FieldTimeSeries(test_name, arch; backend)

            Nt = length(netcdf_JRA55_fts.times)
            @test Oceananigans.OutputReaders.time_indices(netcdf_JRA55_fts) == (1, 2, 3, 4)
            f₁ = view(parent(netcdf_JRA55_fts), :, :, 1, 1)
            f₁ = Array(f₁)

            netcdf_JRA55_fts.backend = Oceananigans.OutputReaders.new_backend(netcdf_JRA55_fts.backend, Nt-2, Nb)
            @test Oceananigans.OutputReaders.time_indices(netcdf_JRA55_fts) == (Nt-2, Nt-1, Nt, 1)
            set!(netcdf_JRA55_fts)

            f₁′ = view(parent(netcdf_JRA55_fts), :, :, 1, 4)
            f₁′ = Array(f₁′)
            @test f₁ == f₁′
        end

        @info "Testing interpolate_field_time_series! on $A..."

        name  = :downwelling_shortwave_radiation
        dates = ClimaOcean.DataWrangling.all_dates(JRA55.RepeatYearJRA55(), name)
        end_date = dates[3]
        JRA55_fts = JRA55FieldTimeSeries(name, arch; end_date)

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

        times = JRA55_fts.times
        boundary_conditions = JRA55_fts.boundary_conditions
        target_fts = FieldTimeSeries{Center, Center, Nothing}(target_grid, times; boundary_conditions)
        interpolate!(target_fts, JRA55_fts)

        # Random regression test
        CUDA.@allowscalar begin
            @test Float32(target_fts[1, 1, 1, 1]) ≈ Float32(222.24313354492188)

            # Only include this if we are filling halo regions within
            # interpolate_field_time_series
            @test Float32(target_fts[Nx + 1, 1, 1, 1]) ≈ Float32(222.24313354492188)
        end

        @test target_fts.times == JRA55_fts.times

        # What else might we test?

        @info "Testing save_field_time_series! on $A..."
        filepath = "JRA55_downwelling_shortwave_radiation_test_$(string(typeof(arch))).jld2" # different filename for each arch so that the CPU and GPU tests do not crash
        ClimaOcean.DataWrangling.save_field_time_series!(target_fts, path=filepath, name="Qsw",
                                                         overwrite_existing = true)
        @test isfile(filepath)

        # Test that we can load the data back
        Qswt = FieldTimeSeries(filepath, "Qsw")
        @test on_architecture(CPU(), parent(Qswt.data)) == on_architecture(CPU(), parent(target_fts.data))
        @test Qswt.times == target_fts.times
        rm(filepath)

        #####
        ##### JRA55 prescribed atmosphere
        #####

        backend    = JRA55NetCDFBackend(2)
        atmosphere = JRA55PrescribedAtmosphere(arch; backend, include_rivers_and_icebergs=false)
        @test atmosphere isa PrescribedAtmosphere
        @test isnothing(atmosphere.auxiliary_freshwater_flux)

        # Test that rivers and icebergs are included in the JRA55 data with the correct frequency
        atmosphere = JRA55PrescribedAtmosphere(arch; backend, include_rivers_and_icebergs=true)
        @test haskey(atmosphere.auxiliary_freshwater_flux, :rivers)
        @test haskey(atmosphere.auxiliary_freshwater_flux, :icebergs)

        rivers_times = atmosphere.auxiliary_freshwater_flux.rivers.times
        pressure_times = atmosphere.pressure.times
        @test rivers_times != pressure_times
        @test length(rivers_times) != length(pressure_times)
        @test rivers_times[2] - rivers_times[1] == 86400

        @info "Testing multi year JRA55 data on $A..."
        dataset = JRA55.MultiYearJRA55()
        dates = ClimaOcean.DataWrangling.all_dates(dataset, :temperature)

        # Test that when date range spans two years both netCDF files are downloaded
        # and concatenated when reading the data.
        start_date = DateTime("1959-01-01T00:00:00") - 15 * Day(1) # sometime in 1958
        end_date   = DateTime("1959-01-01T00:00:00") + 85 * Day(1) # sometime in 1959

        backend = JRA55NetCDFBackend(10)
        Ta = JRA55FieldTimeSeries(:temperature; dataset, start_date, end_date, backend)

        @test Second(end_date - start_date).value ≈ Ta.times[end-1] - Ta.times[1]

        # Test we can access all the data
        for t in eachindex(Ta.times)
            @test Ta[t] isa Field
        end
    end
end
