using Oceananigans
using ClimaOcean
using ClimaOcean.Bathymetry: remove_minor_basins!

@testset "Availability of Bathymetry" begin
    @info "Testing Bathymetry utils..."
    for arch in test_architectures
        grid = LatitudeLongitudeGrid(arch;
                                     size = (100, 100, 10), 
                                     longitude = (0, 100), 
                                     latitude = (0, 50),
                                     z = (-6000, 0))

        # Test that remove_minor_basins!(Z, Inf) does nothing
        bottom_height = regrid_bathymetry(grid)        
        control_bottom_height = deepcopy(bottom_height)
        remove_minor_basins!(bottom_height, Inf)
        @test interior(bottom_height) .== interior(control_bottom_height)

        # Test that remove_minor_basins!(Z, 2) remove the correct number of Basins
        bottom_height = Field{Center, Center, Nothing}(grid)
        control_bottom_height = Field{Center, Center, Nothing}(grid)
        
        # A two basins bathymetry
        bottom(x, y) = - 1000 * Int((x < 10) | (x > 45)) 
        
        set!(bottom_height, bottom)
        set!(control_bottom_height, bottom)

        # This should have not changed anything
        remove_minor_basins!(bottom_height, 2)

        @test interior(bottom_height) .== interior(control_bottom_height)

        # This should have removed the right basin
        remove_minor_basins!(bottom_height, 1)
        
        # The remaning bottom cells that are not immersed should be
        # only (Ny * 20), with a bottom height of 1000, leading to a
        # sum of - Ny * 20 * minimum(bottom_height)
        @test sum(bottom_height) == - grid.Ny * 20 * minimum(bottom_height)

        grid = LatitudeLongitudeGrid(arch;
                                     size = (200, 200, 10), 
                                     longitude = (0, 2), 
                                     latitude = (-10, 50),
                                     z = (-6000, 0))

        control_bottom_height = regrid_bathymetry(grid)
        interpolated_bottom_height = regrid_bathymetry(grid; interpolation_passes = 100)

        # Testing that multiple passes do not change the solution when refining the grid
        @test interior(control_bottom_height) .== interior(interpolated_bottom_height)
    end
end 