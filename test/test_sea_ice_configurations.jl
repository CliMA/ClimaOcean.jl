using ClimaOcean
using ClimaOcean.SeaIceConfigurations
using Oceananigans
using Test

@testset "SeaIceConfigurations method signatures" begin
    @info "Testing SeaIceConfigurations method signatures..."

    # All sea ice constructors accept a positional ocean argument
    for fn in (latitude_longitude_sea_ice,
               half_degree_tripolar_sea_ice,
               one_degree_tripolar_sea_ice,
               sixth_degree_tripolar_sea_ice,
               orca_sea_ice)
        @test hasmethod(fn, Tuple{Any})
    end
end
