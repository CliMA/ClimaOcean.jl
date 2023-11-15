include("runtests_setup.jl")

using ClimaOcean: ECCO2

@testset "ECCO2 fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ecco2_field on $A..."

        temperature_filename = ECCO2.ecco2_file_names[:temperature]
        ecco2_temperature    = ECCO2.ecco2_field(:temperature; architecture=arch)

        @test isfile(temperature_filename)
        rm(temperature_filename)

        @test ecco2_temperature isa Field
        @test ecco2_temperature.grid isa LatitudeLongitudeGrid
        @test topology(ecco2_temperature.grid) == (Periodic, Bounded, Bounded)

        Nx, Ny, Nz = size(ecco2_temperature)
        @test Nx == 1440
        @test Ny == 720
        @test Nz == 50

        ice_thickness_filename = ECCO2.ecco2_file_names[:effective_ice_thickness]
        ecco2_ice_thickness    = ECCO2.ecco2_field(:effective_ice_thickness; architecture=arch)

        @test isfile(ice_thickness_filename)
        rm(ice_thickness_filename)

        @test ecco2_ice_thickness isa Field
        @test ecco2_ice_thickness.grid isa LatitudeLongitudeGrid
        @test topology(ecco2_ice_thickness.grid) == (Periodic, Bounded, Flat)

        Nx, Ny, Nz = size(ecco2_ice_thickness)
        @test Nx == 1440
        @test Ny == 720
        @test Nz == 1
    end
end