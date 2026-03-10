include("runtests_setup.jl")

using Oceananigans
using Statistics
using ClimaOcean

using ClimaOcean.Bathymetry: remove_minor_basins!
using ClimaOcean.Bathymetry: modify_bathymetry_depth!
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


        grid = TripolarGrid(arch; 
                            size = (100, 100, 10),
                            z = (-6000, 0),
                            fold_topology = Oceananigans.Grids.RightFaceFolded)

        @test begin
            regrid_bathymetry(grid)
            true
        end
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

        mask = on_architecture(CPU(), atlantic.mask)

        # Test that the mask is properly bounded
        # Atlantic mask should not include cells in the Pacific
        # (seed point at -170°, 0° should be 0)
        pacific_point_i = findfirst(i -> -175 < i < -165, range(-180, 180, length=360))
        equator_j = 90  # equator for 180 latitude points
        if !isnothing(pacific_point_i)
            @test mask[pacific_point_i, equator_j, 1] == 0
        end
    end
end

# %%
# using KernelAbstractions: get_backend, synchronize

# Assuming the function is available via ClimaOcean or included directly
# using ClimaOcean.Bathymetry: modify_bathymetry_depth!

@testset "Manual Bathymetry Carving (modify_bathymetry_depth!)" begin
    @info "Testing manual bathymetry carving on a global 1-degree grid..."
    
    # Fallback just in case test_architectures isn't globally defined in the test suite
    archs = @isdefined(test_architectures) ? test_architectures : [CPU()]
    
    for arch in archs
        # 1. Setup a global 1-degree grid to precisely match the simulation environment
        #    Range: -180 to 180 degrees (360 cells), -90 to 90 degrees (180 cells)
        grid = LatitudeLongitudeGrid(arch;
                                     size = (360, 180, 1),
                                     longitude = (-180, 180),
                                     latitude = (-90, 90),
                                     z = (0, 1))

        # Initialize flat bathymetry at -1000m
        h = Field{Center, Center, Nothing}(grid)
        set!(h, -1000.0)

        # 2. Prepare 2D Coordinate Arrays using robust mathematical bounds
        #    This completely bypasses any xnode/ynode API quirks across Oceananigans versions
        #    and ensures the coordinate matrices perfectly map to the 1-degree grid.
        lons_2d = zeros(grid.Nx, grid.Ny)
        lats_2d = zeros(grid.Nx, grid.Ny)
        
        for j in 1:grid.Ny, i in 1:grid.Nx
            lons_2d[i, j] = -180.0 + (i - 0.5) * (360.0 / grid.Nx)
            lats_2d[i, j] =  -90.0 + (j - 0.5) * (180.0 / grid.Ny)
        end

        # Move coordinates to the testing architecture (GPU/CPU)
        lons_dev = on_architecture(arch, lons_2d)
        lats_dev = on_architecture(arch, lats_2d)

        # Target the exact center of the global map (near Equator / Prime Meridian)
        target_i, target_j = 180, 90
        target_lon = lons_2d[target_i, target_j]
        target_lat = lats_2d[target_i, target_j]

        # --- TEST 1: Basic Carving & Depth Accuracy ---
        modify_bathymetry_depth!(grid, h, lons_dev, lats_dev, target_lon, target_lat, -100.0; 
                                 rx_km=150.0, ry_km=150.0)
        
        # Pulling Array safely fetches the data back to CPU
        h_cpu = Array(interior(h, :, :, 1))
        
        # Because we targeted the exact cell center, distance = 0, so weight = 1.0
        center_val = h_cpu[target_i, target_j] 
        @test center_val ≈ -100.0 atol=1.0

        # Points far away should remain un-carved
        corner_val = h_cpu[1, 1]
        @test corner_val ≈ -1000.0

        # --- TEST 2: Ellipsoidal Anisotropy (Stretching) ---
        set!(h, -1000.0) # Reset Bathymetry

        # Carve an ellipsoid: Wide in X (400km), Narrow in Y (50km)
        modify_bathymetry_depth!(grid, h, lons_dev, lats_dev, target_lon, target_lat, -100.0; 
                                 rx_km=400.0, ry_km=50.0)

        h_cpu = Array(interior(h, :, :, 1))
        
        # Point A: East (+2 cells at 1.0 deg/cell = 2.0 deg = ~222km). 
        # Inside the 400km rx_km radius -> Should be modified
        val_east = h_cpu[target_i + 2, target_j] 
        @test val_east > -1000.0 

        # Point B: North (+1 cell at 1.0 deg/cell = 1.0 deg = ~111km). 
        # Outside the 50km ry_km radius -> Should remain unmodified
        val_north = h_cpu[target_i, target_j + 1]
        @test val_north ≈ -1000.0

        # --- TEST 3: Smooth Tapering ---
        # Check a point halfway inside the East-West radius (+1 cell = 1.0 deg = ~111km)
        val_mid = h_cpu[target_i + 1, target_j] 
        
        # It should be mathematically blended smoothly between -100 and -1000
        @test -1000.0 < val_mid < -100.0
    end
end