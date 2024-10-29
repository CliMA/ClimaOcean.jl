include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: ECCO_field, metadata_path
using Oceananigans.Grids: topology

using CFTime
using Dates

@testset "ECCO fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ECCO_field on $A..."

        start_date = DateTimeProlepticGregorian(1993, 1, 1)
        end_date = DateTimeProlepticGregorian(1993, 4, 1)
        dates = start_date : Month(1) : end_date

        temperature = ECCOMetadata(:temperature, dates)
        t_restoring = ECCO_restoring_forcing(temperature; timescale = 1000.0)

        ECCO_fts = t_restoring.func.ECCO_fts

        for metadatum in temperature
            @test isfile(metadata_path(metadatum))
        end

        @test ECCO_fts isa FieldTimeSeries
        @test ECCO_fts.grid isa LatitudeLongitudeGrid
        @test topology(ECCO_fts.grid) == (Periodic, Bounded, Bounded)

        Nx, Ny, Nz = size(interior(ECCO_fts))
        Nt = length(ECCO_fts.times)

        @test Nx == size(temperature)[1]
        @test Ny == size(temperature)[2]
        @test Nz == size(temperature)[3]
        @test Nt == size(temperature)[4]
    end
end

@testset "Setting a field with ECCO" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(size=(10, 10, 10), latitude=(-60, -40), longitude=(10, 15), z=(-200, 0))
        field = CenterField(grid)
        set!(field, ECCOMetadata(:temperature)) 
        set!(field, ECCOMetadata(:salinity))
    end
end

@testset "Setting temperature and salinity to ECCO" begin
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(size=(10, 10, 10), latitude=(-60, -40), longitude=(10, 15), z=(-200, 0), halo = (7, 7, 7))
        ocean = ocean_simulation(grid)
        date = DateTimeProlepticGregorian(1993, 1, 1)
        set!(ocean.model, T=ECCOMetadata(:temperature, date), S=ECCOMetadata(:salinity, date))
    end
end
