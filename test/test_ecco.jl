include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: ecco_field, metadata_filename
using Oceananigans.Grids: topology

using CFTime
using Dates

# Make data directory if it doesn't exist
try
    mkdir(ECCO_data_directory)
catch
end

@testset "ECCO fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ecco_field on $A..."

        start_date = DateTimeProlepticGregorian(1993, 1, 1)
        end_date = DateTimeProlepticGregorian(1993, 4, 1)
        dates = start_date : Month(1) : end_date

        temperature = ECCOMetadata(:temperature, dates; dir=ECCO_data_directory)
        t_restoring = ECCO_restoring_forcing(temperature; timescale = 1000.0)

        ecco_fts = t_restoring.func.ecco_fts

        for metadata in temperature
            temperature_filename = metadata_filename(metadata)
            @test isfile(joinpath(ECCO_data_directory, temperature_filename))
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

@testset "setting a field with ECCO" begin
    for arch in test_architectures
        grid  = LatitudeLongitudeGrid(size = (10, 10, 10), latitude = (-60, -40), longitude = (10, 15), z = (-200, 0))
        field = CenterField(grid)
        set!(field, ECCOMetadata(:temperature; dir=ECCO_data_directory)) 
        set!(field, ECCOMetadata(:salinity; dir=ECCO_data_directory))
    end 
end

rm(ECCO_data_directory; recursive=true)