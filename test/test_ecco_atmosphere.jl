include("runtests_setup.jl")

using ClimaOcean.Atmospheres: PrescribedAtmosphere, TwoBandDownwellingRadiation
using ClimaOcean.ECCO: ECCOPrescribedAtmosphere, ECCO4Monthly

@testset "ECCO Prescribed Atmosphere" begin
    for arch in test_architectures
        A = typeof(arch)
        @info "Testing ECCOPrescribedAtmosphere on $A..."

        dataset = ECCO4Monthly()
        start_date = DateTime(1992, 1, 1)
        end_date = DateTime(1992, 3, 1)

        atmosphere = ECCOPrescribedAtmosphere(arch;
                                              dataset,
                                              start_date,
                                              end_date,
                                              time_indices_in_memory = 2)

        @test atmosphere isa PrescribedAtmosphere

        # Test that all expected fields are present
        @test haskey(atmosphere.velocities, :u)
        @test haskey(atmosphere.velocities, :v)
        @test haskey(atmosphere.tracers, :T)
        @test haskey(atmosphere.tracers, :q)
        @test !isnothing(atmosphere.pressure)
        @test !isnothing(atmosphere.downwelling_radiation)
        @test haskey(atmosphere.freshwater_flux, :rain)

        # Test downwelling radiation components
        Qs = atmosphere.downwelling_radiation.shortwave
        Ql = atmosphere.downwelling_radiation.longwave

        @test Qs isa FieldTimeSeries
        @test Ql isa FieldTimeSeries

        # Test that downwelling radiation has the correct sign convention
        # Downwelling radiation should be positive (energy coming down to the surface)
        # This is consistent with JRA55 convention
        CUDA.@allowscalar begin
            # Get some sample values (not all zeros)
            Qs_data = parent(Qs)
            Ql_data = parent(Ql)

            # Longwave radiation should be positive (always some downwelling longwave)
            @test all(Ql_data .>= 0)

            # Typical ranges for radiation (sanity checks)
            # Shortwave: 0 to ~1400 W/m² (solar constant at TOA)
            # Longwave: ~100 to ~500 W/m² (typical atmospheric emission)
            @test maximum(Qs_data) < 1500  # W/m²
            @test maximum(Ql_data) < 600   # W/m²
            @test maximum(Ql_data) > 50    # W/m² - should have some reasonable values
        end

        # Test grid consistency
        @test atmosphere.velocities.u.grid isa LatitudeLongitudeGrid
        @test atmosphere.velocities.u.grid == atmosphere.velocities.v.grid
        @test atmosphere.velocities.u.grid == atmosphere.tracers.T.grid
    end
end
