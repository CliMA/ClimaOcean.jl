include("runtests_setup.jl")

@testset "ECCO2 fields utilities" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ecco2_field on $A..."

        # This should download a file called "RYF.rsds.1990_1991.nc"
        ecco2_temperature = ClimaOcean.ECCO2.ecco2_field(:temperature; architecture=arch)

        @test isfile(ClimaOcean.ECCO2.ecco2_file_names[:temperature])
        rm(ClimaOcean.ECCO2.ecco2_file_names[:temperature])

        @test ecco2_temperature isa Field
        @test ecco2_temperature.grid isa LatitudeLongitudeGrid

        Nx, Ny, Nz = size(ecco2_temperature)
        @test Nx == 1440
        @test Ny == 720
        @test Nz == 50


        # This should download a file called "RYF.rsds.1990_1991.nc"
        ecco2_ice_thickness = ClimaOcean.ECCO2.ecco2_field(:effective_ice_thickness; architecture=arch)

        @test isfile(ClimaOcean.ECCO2.ecco2_file_names[:effective_ice_thickness])
        rm(ClimaOcean.ECCO2.ecco2_file_names[:effective_ice_thickness])

        @test effective_ice_thickness isa Field
        @test effective_ice_thickness.grid isa LatitudeLongitudeGrid
        @test topology(effective_ice_thickness.grid) == (Periodic, Bounded, Flat)

        Nx, Ny, Nz = size(effective_ice_thickness)
        @test Nx == 1440
        @test Ny == 720
        @test Nz == 1
    end
end