include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO4
using Oceananigans.Grids: topology

@testset "ECCO4 fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ecco4_field on $A..."

        temperature_filename = ECCO4.ecco4_file_names[:temperature]
        ecco4_temperature    = ECCO4.ecco4_field(:temperature; architecture=arch)

        @test isfile(temperature_filename)
        rm(temperature_filename)

        @test ecco4_temperature isa Field
        @test ecco4_temperature.grid isa LatitudeLongitudeGrid
        @test topology(ecco4_temperature.grid) == (Periodic, Bounded, Bounded)

        Nx, Ny, Nz = size(ecco4_temperature)
        @test Nx == 1440
        @test Ny == 720
        @test Nz == 50

        ice_thickness_filename = ECCO4.ecco4_file_names[:sea_ice_thickness]
        ecco4_ice_thickness    = ECCO4.ecco4_field(:sea_ice_thickness; architecture=arch)

        @test isfile(ice_thickness_filename)
        rm(ice_thickness_filename)

        @test ecco4_ice_thickness isa Field
        @test ecco4_ice_thickness.grid isa LatitudeLongitudeGrid
        @test topology(ecco4_ice_thickness.grid) == (Periodic, Bounded, Flat)

        Nx, Ny, Nz = size(ecco4_ice_thickness)
        @test Nx == 1440
        @test Ny == 720
        @test Nz == 1
    end
end

@testset "setting a field with ECCO4" begin
    for arch in test_architectures
        grid  = LatitudeLongitudeGrid(size = (10, 10, 10), latitude = (-60, -40), longitude = (-5, 5), z = (-200, 0))
        field = CenterField(grid)
        set!(field, ECCOMetadata(:temperature)) 
        set!(field, ECCOMetadata(:salinity))
    end 
end
