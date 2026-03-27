using ClimaOcean
using ClimaOcean.OceanConfigurations
using ClimaOcean.SeaIceConfigurations
using ClimaOcean.Diagnostics
using Test

@testset "Module exports and reexports" begin
    @info "Testing module exports and reexports..."

    # Core ClimaOcean exports
    @test isdefined(ClimaOcean, :one_degree_tripolar_ocean)
    @test isdefined(ClimaOcean, :half_degree_tripolar_ocean)
    @test isdefined(ClimaOcean, :latitude_longitude_ocean)
    @test isdefined(ClimaOcean, :sixth_degree_tripolar_ocean)
    @test isdefined(ClimaOcean, :Progress)

    # NumericalEarth reexports
    @test isdefined(ClimaOcean, :regrid_bathymetry)
    @test isdefined(ClimaOcean, :ocean_simulation)
    @test isdefined(ClimaOcean, :sea_ice_simulation)

    # Submodule reexports
    @test isdefined(ClimaOcean, :ECCO)
    @test isdefined(ClimaOcean, :ETOPO)
    @test isdefined(ClimaOcean, :JRA55)
    @test isdefined(ClimaOcean, :EN4)
    @test isdefined(ClimaOcean, :GLORYS)

    # Diagnostics exports
    @test isdefined(ClimaOcean.Diagnostics, :MixedLayerDepthField)
    @test isdefined(ClimaOcean.Diagnostics, :compute_report_fields)

    # OceanConfigurations exports
    @test isdefined(OceanConfigurations, :latitude_longitude_ocean)
    @test isdefined(OceanConfigurations, :half_degree_tripolar_ocean)
    @test isdefined(OceanConfigurations, :one_degree_tripolar_ocean)
    @test isdefined(OceanConfigurations, :sixth_degree_tripolar_ocean)
    @test isdefined(OceanConfigurations, :orca_ocean)

    # SeaIceConfigurations exports
    @test isdefined(SeaIceConfigurations, :latitude_longitude_sea_ice)
    @test isdefined(SeaIceConfigurations, :half_degree_tripolar_sea_ice)
    @test isdefined(SeaIceConfigurations, :one_degree_tripolar_sea_ice)
    @test isdefined(SeaIceConfigurations, :sixth_degree_tripolar_sea_ice)
    @test isdefined(SeaIceConfigurations, :orca_sea_ice)
end

@testset "Progress callback" begin
    @info "Testing Progress callback..."

    p = Progress()
    @test p isa Progress
    @test p isa Function
    @test p.wall_time isa Ref{UInt64}
    @test p.wall_time[] > 0

    # Test that Progress(Ref(time_ns())) also works
    t = Ref(time_ns())
    p2 = Progress(t)
    @test p2.wall_time === t
end
