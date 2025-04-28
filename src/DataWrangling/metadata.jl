using CFTime
using Dates
using Base: @propagate_inbounds

struct Metadata{V, D}
    name  :: Symbol
    dataset :: V
    dates :: D
    dir :: String
end

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
- `dataset`: Supported datasets are `ECCO2Monthly()`, `ECCO2Daily()`, `ECCO4Monthly()`, `EN4Monthly(),
             `RepeatYearJRA55()`, or `MultiYearJRA55()`.
- `dates`: The dates of the dataset (`Dates.AbstractDateTime` or `CFTime.AbstractCFDateTime`).
           Note this can either be a range or a vector of dates, representing a time-series.
           For a single date, use [`Metadatum`](@ref).
- `start_date`: If `dates = nothing`, we can prescribe the first date of metadata as a date
                (`Dates.AbstractDateTime` or `CFTime.AbstractCFDateTime`). `start_date` should lie
                within the date range of the dataset. Default: nothing.
- `end_date`: If `dates = nothing`, we can prescribe the last date of metadata as a date
              (`Dates.AbstractDateTime` or `CFTime.AbstractCFDateTime`). `end_date` should lie
                within the date range of the dataset. Default: nothing.
- `dir`: The directory where the dataset is stored.
"""
function Metadata(variable_name;
                  dataset,
                  dates = all_dates(dataset, variable_name),
                  dir = default_download_directory(dataset),
                  start_date = nothing,
                  end_date = nothing)

    if !isnothing(start_date) && !isnothing(end_date)
        dates = compute_native_date_range(dates, start_date, end_date)
    end

    return Metadata(variable_name, dataset, dates, dir)
end

const AnyDateTime  = Union{AbstractCFDateTime, Dates.AbstractDateTime}
const Metadatum{V} = Metadata{V, <:AnyDateTime} where V

"""
    Metadatum(variable_name;
              dataset,
              date = first_date(dataset, variable_name),
              dir = default_download_directory(dataset))

A specialized constructor for a [`Metadata`](@ref) object with a single date, representative of a snapshot in time.
"""
function Metadatum(variable_name;
                   dataset,
                   date = first_date(dataset, variable_name),
                   dir = default_download_directory(dataset))

    date isa Union{CFTime.AbstractCFDateTime, Dates.AbstractDateTime} ||
        throw(ArgumentError("date must be Union{Dates.AbstractDateTime, CFTime.AbstractCFDateTime}"))

    return Metadata(variable_name, dataset, date, dir)
end

# Just the current directory
default_download_directory(dataset) = pwd()

# Default download function for a metadata object, to be extended by each dataset
download_dataset(metadata) = nothing

Base.show(io::IO, metadata::Metadata) =
    print(io, "Metadata:", '\n',
    "├── name: $(metadata.name)", '\n',
    "├── dataset: $(metadata.dataset)", '\n',
    "├── dates: $(metadata.dates)", '\n',
    "└── dir: $(metadata.dir)")

# Treat Metadata as an array to allow iteration over the dates.
Base.length(metadata::Metadata) = length(metadata.dates)
Base.eltype(metadata::Metadata) = Base.eltype(metadata.dates)

# If only one date, it's a single element array
Base.length(metadata::Metadatum) = 1

@propagate_inbounds Base.getindex(m::Metadata, i::Int) = Metadata(m.name, m.dataset, m.dates[i],   m.dir)
@propagate_inbounds Base.first(m::Metadata)            = Metadata(m.name, m.dataset, m.dates[1],   m.dir)
@propagate_inbounds Base.last(m::Metadata)             = Metadata(m.name, m.dataset, m.dates[end], m.dir)

@inline function Base.iterate(m::Metadata, i=1)
    if (i % UInt) - 1 < length(m)
        return Metadata(m.name, m.dataset, m.dates[i], m.dir), i + 1
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

function short_name end

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
#### Some utilities
####

"""
    all_dates(metadata)

Extract all the dates of the given `metadata` formatted using the `DateTime` type.
Note: This methods needs to be extended for any new dataset.
"""
all_dates(metadata) = all_dates(metadata.dataset, metadata.name)

"""
    first_date(dataset)

Extract the first date of the given dataset using the `DateTime` type.
"""
first_date(dataset) = first(all_dates(dataset))

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

"""
    vertical_interfaces(metadata)

Return an array with the vertical interfaces (``z``-faces) of the dataset
that `metadata` corresponds to.
"""
vertical_interfaces(metadata::Metadata{V}) where V =
    error("vertical_interfaces not implemented for $V")

variable_is_three_dimensional(metadata::Metadata{V}) where V =
    error("variable_is_three_dimensional not implemented for $V")

function dataset_latitude_extent end

"""
    empty_field(metadata::Metadata;
                architecture = CPU(),
                horizontal_halo = (7, 7))

Return an empty field of `metadata` on `architecture` and with `horizontal_halo`s.
"""
function empty_field(metadata::Metadata;
                     architecture = CPU(),
                     horizontal_halo = (7, 7))

    Nx, Ny, Nz, _ = size(metadata)
    loc = location(metadata)
    longitude = (0, 360)
    latitude = dataset_latitude_extent(metadata)
    TX, TY = (Periodic, Bounded)

    if variable_is_three_dimensional(metadata)
        TZ = Bounded
        LZ = Center
        z = vertical_interfaces(metadata)
        halo = (horizontal_halo..., 3)
        sz = (Nx, Ny, Nz)
    else # the variable is two-dimensional
        TZ = Flat
        LZ = Nothing
        z = nothing
        halo = horizontal_halo
        sz = (Nx, Ny)
    end

    grid = LatitudeLongitudeGrid(architecture, Float32; halo, longitude, latitude, z,
                                 size = sz,
                                 topology = (TX, TY, TZ))
    return Field{loc...}(grid)
end

struct Celsius end
struct Kelvin end

function dataset_temperature_units end
