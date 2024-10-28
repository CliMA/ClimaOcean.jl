include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: ECCO_field, metadata_path
using ClimaOcean.DataWrangling: NearestNeighborInpainting
using Oceananigans.Grids: topology

using CFTime
using Dates

start_date = DateTimeProlepticGregorian(1993, 1, 1)
end_date = DateTimeProlepticGregorian(1993, 2, 1)
dates = start_date : Month(1) : end_date

# Inpaint only the first two cells inside the missing mask
inpainting = NearestNeighborInpainting(2)

@testset "ECCO fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ECCO_field on $A..."

        temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
        t_restoring = ECCORestoring(temperature; rate = 1 / 1000.0, inpainting)

        ECCO_fts = t_restoring.field_time_series

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

        @test ECCO_fts.times[1]   == ECCO_times(temperature)[1]
        @test ECCO_fts.times[end] == ECCO_times(temperature)[end]
    end
end

@testset "LinearlyTaperedPolarMask" begin

    grid = LatitudeLongitudeGrid(size = (100, 100, 10), 
                             latitude = (-75, 75), 
                            longitude = (0, 360), 
                                    z = (-200, 0),
                                 halo = (6, 6, 6))
    
    φ₁ = grid.φᵃᶜᵃ[1]
    φ₂ = grid.φᵃᶜᵃ[21]
    φ₃ = grid.φᵃᶜᵃ[80]
    φ₄ = grid.φᵃᶜᵃ[100]
    z₁ = grid.zᵃᵃᶜ[6]

    mask = LinearlyTaperedPolarMask(northern = (φ₃, φ₄), 
                                    southern = (φ₁, φ₂), 
                                    z        = (z₁, 0))

    t_restoring = ECCORestoring(:temperature; 
                                dates, 
                                mask, 
                                rate = 1 / 1000.0,
                                inpainting)

    fill!(t_restoring.field_time_series[1], 1.0)
    fill!(t_restoring.field_time_series[2], 1.0)

    T = CenterField(grid)
    fields = (; T)
    clock  = Clock(; time = 0)

    @test t_restoring(1, 1,   10, grid, clock, fields) == t_restoring.rate
    @test t_restoring(1, 11,  10, grid, clock, fields) == t_restoring.rate / 2
    @test t_restoring(1, 21,  10, grid, clock, fields) == 0
    @test t_restoring(1, 80,  10, grid, clock, fields) == 0
    @test t_restoring(1, 90,  10, grid, clock, fields) == t_restoring.rate / 2
    @test t_restoring(1, 100, 10, grid, clock, fields) == t_restoring.rate
    @test t_restoring(1, 1,   5,  grid, clock, fields) == 0
    @test t_restoring(1, 10,  5,  grid, clock, fields) == 0
end

@testset "Timestepping with ECCORestoring" begin
    for arch in test_architectures

        grid  = LatitudeLongitudeGrid(arch; size = (10, 10, 10), 
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

        FT = ECCORestoring(arch, :temperature; dates, rate = 1 / 1000.0, inpainting)
        ocean = ocean_simulation(grid; forcing = (; T = FT))

        @test begin
            time_step!(ocean)
            time_step!(ocean)
            true
        end
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
        set!(ocean.model, T=ECCOMetadata(:temperature, date), S=ECCOMetadata(:salinity, date); inpainting)
    end
end
