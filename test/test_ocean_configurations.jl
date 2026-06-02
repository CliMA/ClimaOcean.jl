using ClimaOcean
using ClimaOcean.OceanConfigurations: vertical_coordinate, henyey_diffusivity, νhb
using Oceananigans
using Test

@testset "OceanConfigurations utilities" begin
    @info "Testing OceanConfigurations utilities..."

    @testset "vertical_coordinate" begin
        z = vertical_coordinate()
        @test z isa Oceananigans.Grids.ExponentialDiscretization
    end

    @testset "henyey_diffusivity" begin
        # At the equator (y=0), should return the equatorial floor value
        κ_eq = henyey_diffusivity(0, 0, 0, 0)
        @test κ_eq == 2e-6

        # At ±30° should be ≈ 1.5e-5
        κ_30 = henyey_diffusivity(0, 30, 0, 0)
        @test κ_30 ≈ 3e-5 * abs(sind(30))
        @test κ_30 ≈ 1.5e-5

        # At higher latitudes, diffusivity should increase
        @test henyey_diffusivity(0, 60, 0, 0) > henyey_diffusivity(0, 30, 0, 0)

        # Symmetric in latitude
        @test henyey_diffusivity(0, 45, 0, 0) ≈ henyey_diffusivity(0, -45, 0, 0)
    end

    @testset "biharmonic viscosity kernel (νhb)" begin
        grid = RectilinearGrid(CPU(); size=(4, 4, 1), x=(0, 1), y=(0, 1), z=(0, 1))
        # νhb computes Az² / λ — just check it returns a finite positive number
        val = νhb(2, 2, 1, grid, Center(), Center(), Center(), nothing, nothing, 1.0)
        @test isfinite(val)
        @test val > 0
    end
end

@testset "OceanConfigurations method signatures" begin
    @info "Testing OceanConfigurations method signatures..."

    # All configuration constructors should accept (arch; kwargs...) signature
    for fn in (latitude_longitude_ocean,
               half_degree_tripolar_ocean,
               one_degree_tripolar_ocean,
               sixth_degree_tripolar_ocean)
        @test hasmethod(fn, Tuple{})
        @test hasmethod(fn, Tuple{CPU})
    end

    @test hasmethod(ClimaOcean.OceanConfigurations.orca_ocean, Tuple{})
    @test hasmethod(ClimaOcean.OceanConfigurations.orca_ocean, Tuple{CPU})
end
