using ClimaOcean
using ClimaOcean.OceanConfigurations: default_half_degree_closure,
                                      default_one_degree_closure,
                                      vertical_coordinate
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

    @testset "vertical_coordinate accepts Nz and depth" begin
        z = vertical_coordinate(; Nz=100, depth=5500)
        @test z isa Oceananigans.Grids.ExponentialDiscretization
    end

    @testset "vertical_coordinate default is unchanged" begin
        z = vertical_coordinate()
        @test z isa Oceananigans.Grids.ExponentialDiscretization
    end
end

@testset "OMIPConfigurations method signatures" begin
    @info "Testing OMIPConfigurations API signatures..."

    @test hasmethod(omip_simulation, Tuple{Symbol})
    @test hasmethod(add_omip_diagnostics!, Tuple{Oceananigans.Simulation})
end
