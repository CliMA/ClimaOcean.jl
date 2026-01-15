module ERA5

export ERA5Hourly, ERA5Monthly

using NCDatasets
using Printf
using Scratch

using Oceananigans.Fields: Center
using ClimaOcean.DataWrangling: Metadata, Metadatum, metadata_path
using Dates
using Dates: DateTime, Day, Month, Hour

import Oceananigans.Fields: location

import ClimaOcean.DataWrangling:
    all_dates,
    dataset_variable_name,
    default_download_directory,
    longitude_interfaces,
    latitude_interfaces,
    z_interfaces,
    metadata_filename,
    inpainted_metadata_path,
    available_variables,
    retrieve_data,
    metadata_path

import Base: eltype

download_ERA5_cache::String = ""

function __init__()
    global download_ERA5_cache = @get_scratch!("ERA5")
end

#####
##### ERA5 Datasets
#####

abstract type ERA5Dataset end

default_download_directory(::ERA5Dataset) = download_ERA5_cache

struct ERA5Hourly <: ERA5Dataset end
struct ERA5Monthly <: ERA5Dataset end

dataset_name(::ERA5Hourly) = "ERA5Hourly"
dataset_name(::ERA5Monthly) = "ERA5Monthly"

# ERA5 has 0.25 degree resolution: 1440 x 721 grid points
Base.size(::ERA5Dataset, variable) = (1440, 721)

# ERA5 reanalysis data available from 1940 to present (we use a practical range here)
all_dates(::ERA5Hourly, var) = range(DateTime("1940-01-01"), stop=DateTime("2024-12-31"), step=Hour(1))
all_dates(::ERA5Monthly, var) = range(DateTime("1940-01-01"), stop=DateTime("2024-12-01"), step=Month(1))

const ERA5Metadata{D} = Metadata{<:ERA5Dataset, D}
const ERA5Metadatum = Metadatum{<:ERA5Dataset}

# ERA5 is a spatially 2D dataset (atmospheric surface variables)
is_three_dimensional(::ERA5Metadata) = false

# Size includes Nz=1 for 2D surface data
Base.size(metadata::ERA5Metadata) = (1440, 721, 1, length(metadata.dates))
Base.size(::ERA5Metadatum) = (1440, 721, 1, 1)

# Variable name mappings from ClimaOcean names to ERA5/era5cli variable names
# The era5cli tool uses these short names for downloading
ERA5_dataset_variable_names = Dict(
    :temperature                     => "2m_temperature",
    :dewpoint_temperature            => "2m_dewpoint_temperature", 
    :eastward_velocity               => "10m_u_component_of_wind",
    :northward_velocity              => "10m_v_component_of_wind",
    :surface_pressure                => "surface_pressure",
    :mean_sea_level_pressure         => "mean_sea_level_pressure",
    :total_precipitation             => "total_precipitation",
    :sea_surface_temperature         => "sea_surface_temperature",
    :downwelling_shortwave_radiation => "surface_solar_radiation_downwards",
    :downwelling_longwave_radiation  => "surface_thermal_radiation_downwards",
    :total_cloud_cover               => "total_cloud_cover",
    :evaporation                     => "evaporation",
    :specific_humidity               => "specific_humidity",
)

# Variables available for download
ERA5_variable_names = keys(ERA5_dataset_variable_names)

available_variables(::ERA5Dataset) = ERA5_dataset_variable_names

dataset_variable_name(metadata::ERA5Metadata) = ERA5_dataset_variable_names[metadata.name]

# NetCDF short variable names (what's actually in the downloaded files)
# These differ from the era5cli parameter names above
ERA5_netcdf_variable_names = Dict(
    :temperature                     => "t2m",
    :dewpoint_temperature            => "d2m", 
    :eastward_velocity               => "u10",
    :northward_velocity              => "v10",
    :surface_pressure                => "sp",
    :mean_sea_level_pressure         => "msl",
    :total_precipitation             => "tp",
    :sea_surface_temperature         => "sst",
    :downwelling_shortwave_radiation => "ssrd",
    :downwelling_longwave_radiation  => "strd",
    :total_cloud_cover               => "tcc",
    :evaporation                     => "e",
    :specific_humidity               => "q",
)

netcdf_variable_name(metadata::ERA5Metadata) = ERA5_netcdf_variable_names[metadata.name]

"""
    retrieve_data(metadata::ERA5Metadatum)

Retrieve ERA5 data from NetCDF file according to `metadata`.
ERA5 is 2D surface data, so we return a 2D array with an added singleton z-dimension.
"""
function retrieve_data(metadata::ERA5Metadatum)
    path = metadata_path(metadata)
    name = netcdf_variable_name(metadata)
    
    ds = NCDatasets.Dataset(path)
    
    # ERA5 is 2D + time, we take the first time step
    # Data shape is typically (lon, lat) or (lon, lat, time)
    raw_data = ds[name]
    ndim = ndims(raw_data)
    
    if ndim == 2
        data_2d = raw_data[:, :]
    elseif ndim == 3
        data_2d = raw_data[:, :, 1]
    else
        error("Unexpected ERA5 data dimensions: $ndim")
    end
    
    close(ds)
    
    # Add singleton z-dimension for 3D field compatibility
    # Return as (Nx, Ny, 1)
    return reshape(data_2d, size(data_2d, 1), size(data_2d, 2), 1)
end

#####
##### Metadata filename construction
#####

function date_str(date::DateTime)
    y = Dates.year(date)
    m = lpad(Dates.month(date), 2, '0')
    return "$(y)-$(m)"
end

start_date_str(date::DateTime) = date_str(date)
end_date_str(date::DateTime) = date_str(date)
start_date_str(dates::AbstractVector) = date_str(first(dates))
end_date_str(dates::AbstractVector) = date_str(last(dates))

colon2dash(s::String) = replace(s, ":" => "-")
underscore_spaces(s::String) = replace(s, " " => "_")

function bbox_strs(::Nothing)
    return "_nothing", "_nothing"
end

function bbox_strs(c)
    first = @sprintf("_%.1f", c[1])
    second = @sprintf("_%.1f", c[2])
    return first, second
end

function metadata_prefix(metadata::ERA5Metadata)
    var = ERA5_dataset_variable_names[metadata.name]
    dataset = dataset_name(metadata.dataset)
    start_date = start_date_str(metadata.dates)
    end_date = end_date_str(metadata.dates)
    bbox = metadata.bounding_box

    if !isnothing(bbox)
        w, e = bbox_strs(bbox.longitude)
        s, n = bbox_strs(bbox.latitude)
        suffix = string(w, e, s, n)
    else
        suffix = ""
    end

    prefix = string(var, "_", dataset, "_", start_date, "_", end_date, suffix)
    prefix = colon2dash(prefix)
    prefix = underscore_spaces(prefix)
    return prefix
end

function metadata_filename(metadata::ERA5Metadatum)
    prefix = metadata_prefix(metadata)
    return string(prefix, ".nc")
end

function metadata_filename(metadata::ERA5Metadata)
    return [metadata_filename(metadatum) for metadatum in metadata]
end

function inpainted_metadata_filename(metadata::ERA5Metadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::ERA5Metadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

#####
##### Grid interfaces
#####

location(::ERA5Metadata) = (Center, Center, Center)

# ERA5 global coverage: 0-360 longitude, -90 to 90 latitude at 0.25 degree resolution
longitude_interfaces(::ERA5Metadata) = (0, 360)
latitude_interfaces(::ERA5Metadata) = (-90, 90)

# ERA5 is a 2D surface dataset, so z is a single level at the surface
z_interfaces(::ERA5Metadata) = (0, 1)

# ERA5 data is stored as Float32
eltype(::ERA5Metadata) = Float32

end # module ERA5

