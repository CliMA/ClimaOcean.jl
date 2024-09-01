include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: ecco_field, metadata_filename
using Oceananigans.Grids: topology

using CFTime
using Dates

@testset "Constructing and ECCO restoring" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ecco_field on $A..."

        start_date = DateTimeProlepticGregorian(1993, 1, 1)
        end_date = DateTimeProlepticGregorian(1993, 4, 1)
        dates = start_date : Month(1) : end_date

        temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
        t_restoring = ECCO_restoring_forcing(temperature; architecture = arch, timescale = 1000.0)

        ecco_fts = t_restoring.func.ecco_fts

        for metadata in temperature
            temperature_filename = metadata_filename(metadata)
            @test isfile(temperature_filename)
        end

        @test ecco_fts isa FieldTimeSeries 
        @test ecco_fts.grid isa LatitudeLongitudeGrid
        @test topology(ecco_fts.grid) == (Periodic, Bounded, Bounded)

        Nx, Ny, Nz = size(interior(ecco_fts))
        Nt = length(ecco_fts.times)

        @test Nx == size(temperature)[1]
        @test Ny == size(temperature)[2]
        @test Nz == size(temperature)[3]
        @test Nt == size(temperature)[4]
    end
end

@testset "Updating an ECCO restoring indices" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ecco_field on $A..."

        start_date = DateTimeProlepticGregorian(1993, 1, 1)
        end_date = DateTimeProlepticGregorian(1993, 4, 1)
        dates = start_date : Month(1) : end_date

        salinity   = ECCOMetadata(:salinity, dates, ECCO4Monthly()) 
        s_salinity = ECCO_restoring_forcing(salinity; architecture = arch, timescale = 1000.0)

        # Check that the backend is set correctly
        @test s_salinity.func.ecco_fts.backend == arch

    end
end