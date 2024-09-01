include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.ECCO: ecco_field, metadata_filename
using Oceananigans.Grids: topology

using CFTime
using Dates

@testset "setting a field with ECCO" begin
    for arch in test_architectures
        grid  = LatitudeLongitudeGrid(size = (10, 10, 10), latitude = (-60, -40), longitude = (10, 15), z = (-200, 0))
        field = CenterField(grid)
        set!(field, ECCOMetadata(:temperature)) 
        set!(field, ECCOMetadata(:salinity))
    end 
end
