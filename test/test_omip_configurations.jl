using ClimaOcean
using ClimaOcean.OceanConfigurations: default_half_degree_closure,
                                      default_one_degree_closure,
                                      vertical_coordinate,
                                      henyey_diffusivity
using ClimaOcean.OMIPConfigurations: omip_simulation, add_omip_diagnostics!
using Oceananigans
using Oceananigans.Units
using Test

@testset "OceanConfigurations parameter API" begin
    @info "Testing parameter-based closure construction..."

    @testset "default_half_degree_closure accepts kwargs" begin
        closure = default_half_degree_closure(;
            κ_skew = 300,
            κ_symmetric = 100,
            biharmonic_timescale = 20days,
        )
        @test length(closure) == 4
        # CATKE
        @test closure[1] isa Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivity
        # Eddy closure with custom parameters
        @test closure[2].κ_skew == 300
        @test closure[2].κ_symmetric == 100
    end

    @testset "default_half_degree_closure preserves defaults" begin
        closure = default_half_degree_closure()
        @test closure[2].κ_skew == 500
        @test closure[2].κ_symmetric == 200
    end

    @testset "default_one_degree_closure accepts kwargs" begin
        closure = default_one_degree_closure(;
            κ_skew = 400,
            κ_symmetric = 150,
        )
        @test closure[2].κ_skew == 400
        @test closure[2].κ_symmetric == 150
    end

    @testset "default_one_degree_closure preserves defaults" begin
        closure = default_one_degree_closure()
        @test closure[2].κ_skew == 500
        @test closure[2].κ_symmetric == 200
    end

    @testset "vertical_coordinate accepts Nz and depth" begin
        z100 = vertical_coordinate(; Nz=100, depth=5500)
        z60  = vertical_coordinate()
        @test z100 isa Oceananigans.Grids.ExponentialDiscretization
        @test z60  isa Oceananigans.Grids.ExponentialDiscretization
    end

    @testset "henyey_diffusivity latitude dependence" begin
        κ_eq = henyey_diffusivity(0, 0, 0, 0)
        κ_30 = henyey_diffusivity(0, 30, 0, 0)
        κ_60 = henyey_diffusivity(0, 60, 0, 0)
        @test κ_eq == 2e-6
        @test κ_30 ≈ 1.5e-5
        @test κ_60 > κ_30
        @test henyey_diffusivity(0, 45, 0, 0) ≈ henyey_diffusivity(0, -45, 0, 0)
    end
end

@testset "OMIPConfigurations API" begin
    @info "Testing OMIPConfigurations API..."

    @testset "omip_simulation method signatures" begin
        @test hasmethod(omip_simulation, Tuple{Symbol})
        # Default config argument
        @test hasmethod(omip_simulation, Tuple{})
    end

    @testset "add_omip_diagnostics! method signature" begin
        @test hasmethod(add_omip_diagnostics!, Tuple{Oceananigans.Simulation})
    end
end
