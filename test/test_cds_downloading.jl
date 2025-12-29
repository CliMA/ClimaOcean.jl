include("runtests_setup.jl")

using CopernicusClimateDataStore
using NCDatasets

using ClimaOcean.DataWrangling.ERA5
using ClimaOcean.DataWrangling.ERA5: ERA5Hourly, ERA5Monthly, ERA5_dataset_variable_names
using ClimaOcean.DataWrangling: metadata_path, download_dataset

using Dates

# Test date: Kyoto Protocol ratification date, February 16, 2005
start_date = DateTime(2005, 2, 16, 12)

@testset "ERA5 data downloading and utilities" begin
    @info "Testing ERA5 downloading and NetCDF file verification..."

    dataset = ERA5Hourly()
    
    # Use a small bounding box to reduce download time
    bounding_box = ClimaOcean.DataWrangling.BoundingBox(longitude=(0, 5), latitude=(40, 45))

    @testset "Download ERA5 temperature data" begin
        variable = :temperature
        metadatum = Metadatum(variable; dataset, bounding_box, date=start_date)

        # Clean up any existing file
        filepath = metadata_path(metadatum)
        isfile(filepath) && rm(filepath; force=true)

        # Download the data
        download_dataset(metadatum)
        @test isfile(filepath)

        # Verify the NetCDF file structure
        ds = NCDataset(filepath)
        
        # Check that it has the expected variable (t2m for 2m_temperature)
        @test haskey(ds, "t2m")
        
        # Check that it has coordinate variables
        @test haskey(ds, "longitude")
        @test haskey(ds, "latitude")
        @test haskey(ds, "time") || haskey(ds, "valid_time")
        
        # Check data dimensions
        lon = ds["longitude"][:]
        lat = ds["latitude"][:]
        @test length(lon) > 0
        @test length(lat) > 0
        
        # Check that data is within expected bounds
        @test minimum(lon) >= -1  # Allow some tolerance
        @test maximum(lon) <= 6
        @test minimum(lat) >= 39
        @test maximum(lat) <= 46
        
        # Check that the temperature data exists and is valid
        t2m = ds["t2m"]
        @test ndims(t2m) >= 2
        
        close(ds)
        
        # Clean up
        rm(filepath; force=true)
    end

    @testset "Availability of ERA5 variables" begin
        # Test that we have defined the key ERA5 variables
        @test haskey(ERA5_dataset_variable_names, :temperature)
        @test haskey(ERA5_dataset_variable_names, :eastward_velocity)
        @test haskey(ERA5_dataset_variable_names, :northward_velocity)
        @test haskey(ERA5_dataset_variable_names, :surface_pressure)
        @test haskey(ERA5_dataset_variable_names, :downwelling_shortwave_radiation)
        @test haskey(ERA5_dataset_variable_names, :downwelling_longwave_radiation)
        
        # Verify variable name mappings
        @test ERA5_dataset_variable_names[:temperature] == "2m_temperature"
        @test ERA5_dataset_variable_names[:eastward_velocity] == "10m_u_component_of_wind"
        @test ERA5_dataset_variable_names[:northward_velocity] == "10m_v_component_of_wind"
    end

    @testset "ERA5 metadata properties" begin
        variable = :temperature
        metadatum = Metadatum(variable; dataset, bounding_box, date=start_date)
        
        # Test metadata properties
        @test metadatum.name == :temperature
        @test metadatum.dataset isa ERA5Hourly
        @test metadatum.dates == start_date
        @test metadatum.bounding_box == bounding_box
        
        # Test size (should be global ERA5 size with 1 time step)
        Nx, Ny, Nz, Nt = size(metadatum)
        @test Nx == 1440  # ERA5 longitude points
        @test Ny == 721   # ERA5 latitude points  
        @test Nz == 1     # 2D surface data
        @test Nt == 1     # Single time step
        
        # Test that ERA5 is correctly identified as 2D
        @test ClimaOcean.DataWrangling.ERA5.is_three_dimensional(metadatum) == false
    end

    @testset "ERA5 Monthly dataset" begin
        monthly_dataset = ERA5Monthly()
        @test monthly_dataset isa ERA5Monthly
        
        # Test that all_dates returns a valid range
        dates = ClimaOcean.DataWrangling.all_dates(monthly_dataset, :temperature)
        @test first(dates) == DateTime("1940-01-01")
        @test step(dates) == Month(1)
    end

    for arch in test_architectures
        A = typeof(arch)

        @testset "Field creation from ERA5 on $A" begin
            variable = :temperature
            metadatum = Metadatum(variable; dataset, bounding_box, date=start_date)

            # Download if not present
            filepath = metadata_path(metadatum)
            isfile(filepath) || download_dataset(metadatum)

            # Create a Field from the downloaded data
            ψ = Field(metadatum, arch)
            @test ψ isa Field

            # ERA5 is 2D data, so field should have Nz=1
            Nx, Ny, Nz = size(ψ)
            @test Nz == 1

            # Verify the field has non-zero data (temperature in Kelvin ~250-310K)
            @allowscalar begin
                @test !all(iszero, interior(ψ))
            end

            # Clean up
            rm(filepath; force=true)
            inpainted_path = ClimaOcean.DataWrangling.inpainted_metadata_path(metadatum)
            isfile(inpainted_path) && rm(inpainted_path; force=true)
        end

        @testset "Setting a field from ERA5 metadata on $A" begin
            variable = :temperature
            metadatum = Metadatum(variable; dataset, bounding_box, date=start_date)

            # Download if not present
            filepath = metadata_path(metadatum)
            isfile(filepath) || download_dataset(metadatum)

            # Create a target grid matching the bounding box region
            grid = LatitudeLongitudeGrid(arch;
                                         size = (10, 10, 1),
                                         latitude = (40, 45),
                                         longitude = (0, 5),
                                         z = (0, 1))

            field = CenterField(grid)

            # Set the field from metadata
            set!(field, metadatum)

            # Verify the field was set with non-zero data
            @allowscalar begin
                @test !all(iszero, interior(field))
            end

            # Clean up
            rm(filepath; force=true)
            inpainted_path = ClimaOcean.DataWrangling.inpainted_metadata_path(metadatum)
            isfile(inpainted_path) && rm(inpainted_path; force=true)
        end
    end
end
