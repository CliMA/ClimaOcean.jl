using Oceananigans
using ClimaOcean

@testset "Availability of Bathymetry" begin
    @info "Testing that we can download the Bathyemtry..."
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size = (100, 100, 10), 
                                     longitude = (0, 100), 
                                     latitude = (-10, 50),
                                     z = (-6000, 0))

        # Just a simple test
        bottom_height = regrid_bathymetry(grid)
    end
end 