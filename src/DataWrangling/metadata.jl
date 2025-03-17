using CFTime
using Dates
using Base: @propagate_inbounds

struct Metadata{D, V}
    name  :: Symbol
    dates :: D
    dataset :: V
    dir :: String
end

"""
   Metadata(variable_name;
            dataset,
            dates = all_dates(dataset, variable_name),
            dir = default_download_directory(dataset))

Metadata holding a specific dataset information.

Arguments
=========
- `variable_name`: a symbol representing the name of the variable (for example, `:temperature`,
                   `:salinity`, `:u_velocity`, etc)

Keyword Arguments
=================
- `dataset`: The dataset of the dataset. Supported datasets are `ECCO2Monthly()`, `ECCO2Daily()`,
             `ECCO4Monthly()`, `JRA55RepeatYear()`, or `JRA55MultipleYears()`.
- `dates`: The dates of the dataset, in a `AbstractCFDateTime` format. Note this can either be a range
           or a vector of dates, representing a time-series. For a single date, use [`Metadatum`](@ref).
- `dir`: The directory where the dataset is stored.
"""
function Metadata(variable_name;
                  dataset,
                  dates=all_dates(dataset, variable_name)[1:1],
                  dir=default_download_directory(dataset))

    return Metadata(variable_name, dates, dataset, dir)
end

const AnyDateTime = Union{AbstractCFDateTime, Dates.AbstractDateTime}
const Metadatum   = Metadata{<:AnyDateTime}

"""
    Metadatum(variable_name;
              dataset,
              date=first_date(dataset, variable_name),
              dir=default_download_directory(dataset))

A specialized constructor for a [`Metadata`](@ref) object with a single date, representative of a snapshot in time.
"""
function Metadatum(variable_name;
                   dataset,
                   date=first_date(dataset, variable_name),
                   dir=default_download_directory(dataset))

    # TODO: validate that `date` is actually a single date?
    return Metadata(variable_name, date, dataset, dir)
end

# Just the current directory
default_download_directory(dataset) = pwd()

# Default download function for a metadata object, to be extended by each dataset
download_dataset(metadata) = nothing

Base.show(io::IO, metadata::Metadata) =
    print(io, "ECCOMetadata:", '\n',
    "├── name: $(metadata.name)", '\n',
    "├── dates: $(metadata.dates)", '\n',
    "├── dataset: $(metadata.dataset)", '\n',
    "└── data directory: $(metadata.dir)")

# Treat Metadata as an array to allow iteration over the dates.
Base.length(metadata::Metadata) = length(metadata.dates)
Base.eltype(metadata::Metadata) = Base.eltype(metadata.dates)

# If only one date, it's a single element array
Base.length(metadata::Metadatum) = 1

@propagate_inbounds Base.getindex(m::Metadata, i::Int) = Metadata(m.name, m.dates[i],   m.dataset, m.dir)
@propagate_inbounds Base.first(m::Metadata)            = Metadata(m.name, m.dates[1],   m.dataset, m.dir)
@propagate_inbounds Base.last(m::Metadata)             = Metadata(m.name, m.dates[end], m.dataset, m.dir)

@inline function Base.iterate(m::Metadata, i=1)
    if (i % UInt) - 1 < length(m)
        return Metadata(m.name, m.dates[i], m.dataset, m.dir), i + 1
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

metadata_path(metadata) = joinpath(metadata.dir, metadata_filename(metadata))

"""
    native_times(metadata; start_time=first(metadata).dates)

Extract the time values from the given `metadata` and calculate the time difference
from the `start_time` and return as an array of time differences in seconds.

Arguments
=========
- `metadata`: The metadata containing the date information.
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
#### Some utilities
####

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

"""
    compute_native_date_range(native_dates, start_date, end_date)

Compute the range of dates that fall within the specified start and end date.
"""
function compute_native_date_range(native_dates, start_date, end_date)
    if last(native_dates) < end_date
        @warn "`end_date` ($end_date) is after the last date in the dataset $(last(native_dates))"
    end

    if last(native_dates) < start_date
       throw(ArgumentError("`start_date` ($start_date) is after the last date in the dataset $(last(native_dates))"))
    end

    start_idx = findfirst(x -> x ≥ start_date, native_dates)
    end_idx   = findfirst(x -> x ≥ end_date, native_dates)
    
    start_idx = start_idx > 1 ? start_idx - 1 : start_idx
    end_idx   = isnothing(end_idx) ? length(native_dates) : end_idx

    return native_dates[start_idx:end_idx]
end
