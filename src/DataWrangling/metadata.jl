using CFTime
using Dates
using Base: @propagate_inbounds
using Oceananigans.Utils: prettysummary

struct BoundingBox{X, Y, Z}
    longitude :: X
    latitude :: Y 
    z :: Z
end

latitude_summary(::Nothing) = "latitude=nothing"
longitude_summary(::Nothing) = "longitude=nothing"
latitude_summary(lat) = string("latitude=(", prettysummary(lat[1]), ", ", prettysummary(lat[2]), ")")
longitude_summary(lon) = string("longitude=(", prettysummary(lon[1]), ", ", prettysummary(lon[2]), ")")
Base.summary(bbox::BoundingBox) = string("BoundingBox(",
                                         longitude_summary(bbox.longitude), ", ",
                                         latitude_summary(bbox.latitude), ")")

"""
    BoundingBox(; longitude=nothing, latitude=nothing, z=nothing)

Create a bounding box with `latitude`, `longitude`, and `z` bounds on the sphere.
"""
BoundingBox(; longitude=nothing, latitude=nothing, z=nothing) =
    BoundingBox(longitude, latitude, z)

struct Metadata{V, D, B}
    name  :: Symbol
    dataset :: V
    dates :: D
    bounding_box :: B
    dir :: String
end

is_three_dimensional(::Metadata) = true
z_interfaces(md::Metadata) = z_interfaces(md.dataset)
longitude_interfaces(md::Metadata) = longitude_interfaces(md.dataset)
latitude_interfaces(md::Metadata) = latitude_interfaces(md.dataset)

"""
    Metadata(variable_name;
             dataset,
             dates = all_dates(dataset, variable_name),
             dir = default_download_directory(dataset))

Metadata holding a specific dataset information.

Argument
========
- `variable_name`: a symbol representing the name of the variable (for example, `:temperature`,
                   `:salinity`, `:u_velocity`, etc)

Keyword Arguments
=================

- `dataset`: Supported datasets are `ECCO2Monthly()`, `ECCO2Daily()`, `ECCO4Monthly()`, `EN4Monthly()`,
             `RepeatYearJRA55()`, and `MultiYearJRA55()`.

- `dates`: The dates of the dataset (`Dates.AbstractDateTime` or `CFTime.AbstractCFDateTime`).
           Note that `dates` can either be a range or a vector of dates, representing a time-series.
           For a single date, use [`Metadatum`](@ref).

- `start_date`: If `dates = nothing`, we can prescribe the first date of metadata as a date
                (`Dates.AbstractDateTime` or `CFTime.AbstractCFDateTime`). `start_date` should lie
                within the date range of the dataset. Default: nothing.

- `end_date`: If `dates = nothing`, we can prescribe the last date of metadata as a date
              (`Dates.AbstractDateTime` or `CFTime.AbstractCFDateTime`). `end_date` should lie
                within the date range of the dataset. Default: nothing.

- `bounding_box`: Specifies the bounds of the dataset. See [`BoundingBox`](@ref).

- `dir`: The directory where the dataset is stored.
"""
function Metadata(variable_name;
                  dataset,
                  dates = all_dates(dataset, variable_name),
                  dir = default_download_directory(dataset),
                  bounding_box = nothing,
                  start_date = nothing,
                  end_date = nothing)

    if !isnothing(start_date) && !isnothing(end_date)
        dates = compute_native_date_range(dates, start_date, end_date)
    end

    return Metadata(variable_name, dataset, dates, bounding_box, dir)
end

const AnyDateTime  = Union{AbstractCFDateTime, Dates.AbstractDateTime}
const Metadatum{V} = Metadata{V, <:AnyDateTime} where V

function Base.size(metadata::Metadata)
    Nx, Ny, Nz = size(metadata.dataset, metadata.name)
    if metadata.dates isa AbstractArray
        Nt = length(metadata.dates)
    else
        Nt = 1
    end
    return (Nx, Ny, Nz, Nt)
end

"""
    Metadatum(variable_name;
              dataset,
              bounding_box = nothing,
              date = first_date(dataset, variable_name),
              dir = default_download_directory(dataset))

A specialized constructor for a [`Metadata`](@ref) object with a single date, representative of a snapshot in time.
"""
function Metadatum(variable_name;
                   dataset,
                   bounding_box = nothing,
                   date = first_date(dataset, variable_name),
                   dir = default_download_directory(dataset))

    return Metadata(variable_name, dataset, date, bounding_box, dir)
end

datestr(md::Metadata) = string(first(md.dates), "--", last(md.dates))
datestr(md::Metadatum) = string(md.dates)
datasetstr(md::Metadata) = string(md.dataset)
metaprefix(md::Metadata) = string("Metadata{", md.dataset, "}")

function Base.show(io::IO, metadata::Metadata)
    print(io, "Metadata:", '\n',
    "├── name: $(metadata.name)", '\n',
    "├── dataset: $(metadata.dataset)", '\n',
    "├── dates: $(metadata.dates)", '\n')

    bbox = metadata.bounding_box
    if !isnothing(bbox)
        print(io, "├── bounding_box: ", summary(bbox), '\n')
    end

    print(io, "└── dir: $(metadata.dir)")
end

# Treat Metadata as an array to allow iteration over the dates.
Base.length(metadata::Metadata) = length(metadata.dates)
Base.eltype(metadata::Metadata) = Float32

Base.summary(md::Metadata) = string(metaprefix(md),
                                    "{", datasetstr(md), "} of ",
                                    md.name, " for ", datestr(md))

# If only one date, it's a single element array
Base.length(metadata::Metadatum) = 1

@propagate_inbounds Base.getindex(m::Metadata, i::Int) = Metadata(m.name, m.dataset, m.dates[i],   m.bounding_box, m.dir)
@propagate_inbounds Base.first(m::Metadata)            = Metadata(m.name, m.dataset, m.dates[1],   m.bounding_box, m.dir)
@propagate_inbounds Base.last(m::Metadata)             = Metadata(m.name, m.dataset, m.dates[end], m.bounding_box, m.dir)

@inline function Base.iterate(m::Metadata, i=1)
    if (i % UInt) - 1 < length(m)
        return Metadata(m.name, m.dataset, m.dates[i], m.bounding_box, m.dir), i + 1
    else
        return nothing
    end
end

# Implementation for 1 date
Base.axes(metadata::Metadatum)    = 1
Base.first(metadata::Metadatum)   = metadata
Base.last(metadata::Metadatum)    = metadata
Base.iterate(metadata::Metadatum) = (metadata, nothing)
Base.iterate(::Metadatum, ::Any)  = nothing

metadata_path(metadata::Metadatum) = joinpath(metadata.dir, metadata_filename(metadata))
metadata_path(metadata::Metadata) = [metadata_path(metadatum) for metadatum in metadata]

"""
    native_times(metadata; start_time=first(metadata).dates)

Extract the time values from the given `metadata`, calculate the time difference
from the `start_time`, and return an array of time differences in seconds.

Argument
========
- `metadata`: The metadata containing the date information.

Keyword Argument
================
- `start_time`: The start time for calculating the time difference. Defaults to the first date in the metadata.
"""
function native_times(metadata; start_time=first(metadata).dates)
    times = zeros(length(metadata))
    for (t, data) in enumerate(metadata)
        date = data.dates
        delta = date - start_time
        delta = Second(delta).value
        times[t] = delta
    end

    return times
end

####
#### Metadata interface
####

"""
    download_dataset(metadata)

Download the dataset file(s) specified by `metadata` and return the path(s) to the downloaded file.
"""
function download_dataset end

"""
    default_download_directory(dataset)

Return the default directory to which `dataset` is downloaded.
"""
function default_download_directory end

"""
    dataset_variable_name(metadata)

Return the name used for the variable metadata.name in its raw dataset file.
"""
function dataset_variable_name end

"""
    all_dates(metadata)

Extract all the dates of the given `metadata` formatted using the `DateTime` type.
Note: This methods needs to be extended for any new dataset.
"""
all_dates(metadata) = all_dates(metadata.dataset, metadata.name)

"""
    first_date(dataset, variable_name)

Extracts the first date of the given dataset and variable name formatted using the `DateTime` type.
"""
first_date(dataset, variable_name) = first(all_dates(dataset, variable_name))

"""
    last_date(dataset, variable_name)

Extracts the last date of the given dataset and variable name formatted using the `DateTime` type.
"""
last_date(dataset, variable_name) = last(all_dates(dataset, variable_name))

"""
    metadata_filename(metadata)

File names of metadata containing multiple dates. The specific version for a `Metadatum` object is
extended in the data specific modules.
"""
metadata_filename(metadata) = [metadata_filename(metadatum) for metadatum in metadata]

struct Celsius end
struct Kelvin end

temperature_units(metadata) = Celsius()

#####
##### Utilities
#####

"""
    compute_native_date_range(native_dates, start_date, end_date)

Compute the range of `native_dates` that fall within the specified `start_date` and `end_date`.
"""
function compute_native_date_range(native_dates, start_date, end_date)
    if last(native_dates) < end_date
        @warn "`end_date` ($end_date) is after the last date in the dataset $(last(native_dates))"
    end

    if last(native_dates) < start_date
       throw(ArgumentError("`start_date` ($start_date) is after the last date in the dataset $(last(native_dates))"))
    end

    start_idx = findfirst(x -> x ≥ start_date, native_dates)
    end_idx   = findfirst(x -> x ≥ end_date,   native_dates)
    start_idx = (start_idx > 1 && native_dates[start_idx] > start_date) ? start_idx - 1 : start_idx
    end_idx   = isnothing(end_idx) ? length(native_dates) : end_idx

    return native_dates[start_idx:end_idx]
end
