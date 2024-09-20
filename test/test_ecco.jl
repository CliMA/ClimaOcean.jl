include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: ECCO_field, metadata_filename, ECCO_times
using Oceananigans.Grids: topology

using CFTime
using Dates

start_date = DateTimeProlepticGregorian(1993, 1, 1)
end_date = DateTimeProlepticGregorian(1993, 2, 1)
dates = start_date : Month(1) : end_date

@testset "ECCO fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ECCO_field on $A..."

        temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
        t_restoring = ECCORestoring(temperature; rate = 1 / 1000.0)

        ECCO_fts = t_restoring.ECCO_fts

        for metadata in temperature
            temperature_filename = metadata_filename(metadata)
            @test isfile(temperature_filename)
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

        @test ECCO_fts.times[1]   == ECCO_times(temperature[1])
        @test ECCO_fts.times[end] == ECCO_times(temperature[end])
    end
end

@testset "LatitudinallyTaperedPolarMask" begin

    grid = LatitudeLongitudeGrid(size = (100, 100, 10), latitude = (-75, 75), longitude = (0, 360), z = (-200, 0))
    
    φ₁ = grid.φᵃᶜᵃ[1]
    φ₂ = grid.φᵃᶜᵃ[20]
    φ₃ = grid.φᵃᶜᵃ[80]
    φ₄ = grid.φᵃᶜᵃ[100]
    z₁ = grid.zᵃᵃᶜ[6]

    mask = LatitudinallyTaperedPolarMask(northern_edges = (φ₃, φ₄), 
                                         southern_edges = (φ₁, φ₂), 
                                         z_edges = (z₁, 0))

    t_restoring = ECCORestoring(:temperature, CPU(); dates, mask, rate = 1 / 1000.0)

    fill!(t_restoring.ECCO_fts[1], 1.0)
    fill!(t_restoring.ECCO_fts[2], 1.0)

    T = CenterField(grid)
    fields = (; T)
    clock  = Clock(; time = 0)

    @test t_restoring(1, 1,   10, grid, clock, fields) == t_restoring.rate
    @test t_restoring(1, 20,  10, grid, clock, fields) == 0
    @test t_restoring(1, 80,  10, grid, clock, fields) == 0
    @test t_restoring(1, 100, 10, grid, clock, fields) == t_restoring.rate
    @test t_restoring(1, 1,   5,  grid, clock, fields) == 0
    @test t_restoring(1, 10,  5,  grid, clock, fields) == 0
end

@testset "setting a field with ECCO" begin
    for arch in test_architectures
        grid  = LatitudeLongitudeGrid(arch; size = (10, 10, 10), latitude = (-60, -40), longitude = (10, 15), z = (-200, 0))
        field = CenterField(grid)
        
        @test begin
            set!(field, ECCOMetadata(:temperature)) 
            set!(field, ECCOMetadata(:salinity))
            true
        end

        FT = ECCORestoring(:temperature, arch; rate = 1 / 1000.0)
        ocean = ocean_simulation(grid; forcing = (; T = FT))

        @test begin
            time_step!(ocean)
            time_step!(ocean)
            true
        end
    end 
end
