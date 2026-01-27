include("runtests_setup.jl")

using Oceananigans
using Statistics
using ClimaOcean

using ClimaOcean.Bathymetry: remove_minor_basins!, label_ocean_basins, find_label_at_point
using ClimaOcean.Bathymetry: OceanBasinMask, atlantic_ocean_mask, pacific_ocean_mask
using ClimaOcean.Bathymetry: Barrier, ATLANTIC_OCEAN_BARRIERS
using ClimaOcean.DataWrangling.ETOPO

@testset "Bathymetry construction and smoothing" begin
    @info "Testing Bathymetry construction and smoothing..."
    for arch in test_architectures
        ETOPOmetadata = Metadatum(:bottom_height, dataset=ETOPO2022())

        # Testing downloading
        ClimaOcean.DataWrangling.download_dataset(ETOPOmetadata)
        @test isfile(metadata_path(ETOPOmetadata))

        grid = LatitudeLongitudeGrid(arch;
                                     size = (100, 100, 10),
                                     longitude = (0, 100),
                                     latitude = (0, 50),
                                     z = (-6000, 0))

        # Test that remove_minor_basins!(Z, Inf) does nothing
        control_bottom_height = regrid_bathymetry(grid, ETOPOmetadata)
        bottom_height = deepcopy(control_bottom_height)
        @test_throws ArgumentError remove_minor_basins!(bottom_height, Inf)

        # A fictitiously large number which should presumably keep all the basins
        remove_minor_basins!(bottom_height, 10000000)
        @test parent(bottom_height) == parent(control_bottom_height)

        # Test that remove_minor_basins!(Z, 2) remove the correct number of Basins
        bottom_height = Field{Center, Center, Nothing}(grid)
        control_bottom_height = Field{Center, Center, Nothing}(grid)

        # A two-basins bathymetry
        bottom(x, y) = - 1000 * Int((x < 10) | (x > 50))

        set!(bottom_height, bottom)
        set!(control_bottom_height, bottom)

        # This should have not changed anything
        remove_minor_basins!(bottom_height, 2)
        @test parent(bottom_height) == parent(control_bottom_height)

        # This should have removed the left basin
        remove_minor_basins!(bottom_height, 1)

        # The remaining bottom cells that are not immersed should be only on the right hand side
        # The left half of the domain should be fully immersed, i.e., bottom == 0
        @test sum(view(bottom_height, 1:50, :, 1)) == 0

        # While the right side should be not immersed, with a mean bottom depth
        # of -1000 meters
        @test mean(view(bottom_height, 51:100, :, 1)) == -1000

        grid = LatitudeLongitudeGrid(arch;
                                     size = (200, 200, 10),
                                     longitude = (0, 100),
                                     latitude = (-10, 50),
                                     z = (-6000, 0))

        control_bottom_height = regrid_bathymetry(grid)
        interpolated_bottom_height = regrid_bathymetry(grid; interpolation_passes=10)

        # Testing that multiple passes _do_ change the solution when coarsening the grid
        @test parent(control_bottom_height) != parent(interpolated_bottom_height)
    end
end

@testset "Barrier geometry" begin
    @info "Testing barrier geometry utilities..."

    # Test Barrier construction with explicit bounds
    barrier = Barrier(-10.0, 10.0, -5.0, 5.0)
    @test barrier.west == -10.0
    @test barrier.east == 10.0
    @test barrier.south == -5.0
    @test barrier.north == 5.0

    # Test Barrier with keyword arguments
    barrier_kw = Barrier(west=-10.0, east=10.0, south=-5.0, north=5.0)
    @test barrier_kw.west == -10.0
    @test barrier_kw.east == 10.0

    # Test meridional barrier constructor (3 args + width)
    meridional = Barrier(20.0, -36.0, -30.0)  # longitude, south, north
    @test meridional.west == 19.0   # 20 - 2/2
    @test meridional.east == 21.0   # 20 + 2/2
    @test meridional.south == -36.0
    @test meridional.north == -30.0

    # Test meridional barrier with custom width
    meridional_wide = Barrier(20.0, -36.0, -30.0; width=4.0)
    @test meridional_wide.west == 18.0
    @test meridional_wide.east == 22.0

    # Test LatitudeBand
    band = Barrier(-180.0, 180.0, -60.0, -55.0)
    @test band.west == -180.0
    @test band.east == 180.0
    @test band.south == -60.0
    @test band.north == -55.0
end

@testset "Ocean basin labeling with barriers" begin
    @info "Testing ocean basin labeling with barriers..."

    for arch in test_architectures
        # Create a global grid
        grid = LatitudeLongitudeGrid(arch;
                                     size = (90, 45, 10),
                                     longitude = (-180, 180),
                                     latitude = (-90, 90),
                                     z = (-6000, 0))

        # Regrid real bathymetry onto this grid
        bottom_height = regrid_bathymetry(grid)
        ibg = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

        # Test labeling without barriers - all oceans should have the same label
        # (they're connected via the Southern Ocean)
        labels_no_barrier = label_ocean_basins(ibg)

        # Find labels at Atlantic and Pacific seed points
        atlantic_label = find_label_at_point(labels_no_barrier, ibg, -30.0, 0.0)
        pacific_label = find_label_at_point(labels_no_barrier, ibg, -170.0, 0.0)

        # Without barriers, Atlantic and Pacific should have the same label
        # (connected via Southern Ocean)
        @test atlantic_label == pacific_label
        @test atlantic_label > 0  # Should find a valid basin

        # Test labeling with barriers - oceans should be separated
        labels_with_barrier = label_ocean_basins(ibg; barriers=ATLANTIC_OCEAN_BARRIERS)

        # With barriers, Atlantic should still be found
        atlantic_label_with_barrier = find_label_at_point(labels_with_barrier, ibg, -30.0, 0.0)
        @test atlantic_label_with_barrier > 0
    end
end

@testset "OceanBasinMask creation" begin
    @info "Testing OceanBasinMask creation..."

    for arch in test_architectures
        # Create a global grid at 1° resolution (needed to properly resolve
        # Central America and separate Atlantic from Pacific)
        grid = LatitudeLongitudeGrid(arch;
                                     size = (360, 180, 10),
                                     longitude = (-180, 180),
                                     latitude = (-90, 90),
                                     z = (-6000, 0))

        bottom_height = regrid_bathymetry(grid)
        ibg = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

        # Test atlantic_ocean_mask creation
        atlantic = atlantic_ocean_mask(ibg)
        @test atlantic isa OceanBasinMask
        @test sum(atlantic.mask) > 0  # Should have some ocean cells

        # Test that the mask is properly bounded
        # Atlantic mask should not include cells in the Pacific
        # (seed point at -170°, 0° should be 0)
        pacific_point_i = findfirst(i -> -175 < i < -165, range(-180, 180, length=360))
        equator_j = 90  # equator for 180 latitude points
        if !isnothing(pacific_point_i)
            @test atlantic.mask[pacific_point_i, equator_j, 1] == 0
        end
    end
end
