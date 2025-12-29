module ClimaOceanCopernicusClimateDataStoreExt

using ClimaOcean
using CopernicusClimateDataStore

using Oceananigans
using Oceananigans.DistributedComputations: @root

using Dates
using ClimaOcean.DataWrangling.ERA5: ERA5Metadata, ERA5Metadatum, ERA5_dataset_variable_names

import ClimaOcean.DataWrangling: download_dataset

"""
    download_dataset(metadata::ERA5Metadata; kwargs...)

Download ERA5 data for each date in the metadata, returning paths to downloaded files.
"""
function download_dataset(metadata::ERA5Metadata; kwargs...)
    paths = Array{String}(undef, length(metadata))
    for (m, metadatum) in enumerate(metadata)
        paths[m] = download_dataset(metadatum; kwargs...)
    end
    return paths
end

"""
    download_dataset(meta::ERA5Metadatum; skip_existing=true, kwargs...)

Download ERA5 data for a single date/time using the CopernicusClimateDataStore package.

The download is performed using `era5cli` through the CopernicusClimateDataStore package.

# Keyword Arguments
- `skip_existing`: Skip download if the file already exists (default: `true`).
- Additional keyword arguments are passed to `CopernicusClimateDataStore.hourly`.

# Environment Setup
Before downloading, you must:
1. Create an account at https://cds.climate.copernicus.eu/
2. Accept the Terms of Use for the ERA5 dataset on the dataset page
3. Set up your API credentials in `~/.cdsapirc`

See https://cds.climate.copernicus.eu/how-to-api for details.
"""
function download_dataset(meta::ERA5Metadatum;
                          skip_existing = true,
                          threads = 1,
                          additional_kw...)

    output_directory = meta.dir
    output_filename = ClimaOcean.DataWrangling.metadata_filename(meta)
    output_path = joinpath(output_directory, output_filename)

    # Skip if file already exists
    if skip_existing && isfile(output_path)
        return output_path
    end

    # Ensure output directory exists
    mkpath(output_directory)

    # Get the ERA5 variable name
    variable_name = ERA5_dataset_variable_names[meta.name]

    # Extract date information
    date = meta.dates
    year = Dates.year(date)
    month = Dates.month(date)
    day = Dates.day(date)
    hour = Dates.hour(date)

    # Build area constraint from bounding box
    area = build_era5_area(meta.bounding_box)

    # Build output prefix (filename without extension)
    output_prefix = first(splitext(output_filename))

    # Perform the download using era5cli via CopernicusClimateDataStore
    @root begin
        downloaded_files = CopernicusClimateDataStore.hourly(;
            variables = variable_name,
            startyear = year,
            months = month,
            days = day,
            hours = hour,
            area = area,
            format = "netcdf",
            outputprefix = output_prefix,
            overwrite = !skip_existing,
            threads = threads,
            splitmonths = false,
            directory = output_directory,
            additional_kw...
        )

        # era5cli generates its own filename suffix, so rename to our expected name
        if !isempty(downloaded_files)
            downloaded_file = first(downloaded_files)
            if downloaded_file != output_path && isfile(downloaded_file)
                mv(downloaded_file, output_path; force=true)
            end
        end
    end

    return output_path
end

#####
##### Area/bounding box utilities
#####

build_era5_area(::Nothing) = nothing

const BBOX = ClimaOcean.DataWrangling.BoundingBox

function build_era5_area(bbox::BBOX)
    # ERA5/era5cli uses (lat_max, lon_min, lat_min, lon_max) ordering
    # BoundingBox has longitude = (west, east), latitude = (south, north)
    
    lon = bbox.longitude
    lat = bbox.latitude

    if isnothing(lon) || isnothing(lat)
        return nothing
    end

    lon_min = lon[1]  # west
    lon_max = lon[2]  # east
    lat_min = lat[1]  # south
    lat_max = lat[2]  # north

    # Return in era5cli order: (lat_max, lon_min, lat_min, lon_max)
    return (lat = (lat_min, lat_max), lon = (lon_min, lon_max))
end

end # module ClimaOceanCopernicusClimateDataStoreExt

