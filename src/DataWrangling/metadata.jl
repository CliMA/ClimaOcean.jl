using CFTime
using Dates
using Base: @propagate_inbounds

struct Metadata{D, V} 
    name  :: Symbol
    dates :: D
    version :: V
    dir :: String
end

"""
   Metadata(variable_name;
            version,
            dates = all_dates(name),
            dir = default_download_folder(version))

Metadata holding a specific dataset information.

Arguments
=========
- `variable_name`: a symbol representing the name of the variable (for example :temperature, :salinity, :u_velocity, etc...)

Keyword Arguments
=================
- `dates`: The dates of the dataset, in a `AbstractCFDateTime` format.. Note this can either be a single date,
           representing a snapshot, or a range of dates, representing a time-series.
- `version`: The version of the dataset. Supported versions are `ECCO2Monthly()`, `ECCO2Daily()`, `ECCO4Monthly()`,
             `JRA55RepeatYear()`, or `JRA55MultipleYears()`.
- `dir`: The directory where the dataset is stored.
"""
function Metadata(variable_name;
                  version,
                  dates = all_dates(version),
                  dir = default_download_folder(version))

    return Metadata(variable_name, dates, version, dir)
end

default_download_folder(version) = "./"

Base.show(io::IO, metadata::Metadata) = 
    print(io, "ECCOMetadata:", '\n',
    "├── name: $(metadata.name)", '\n',
    "├── dates: $(metadata.dates)", '\n',
    "├── version: $(metadata.version)", '\n',
    "└── data directory: $(metadata.dir)")

# Treat Metadata as an array to allow iteration over the dates.
Base.length(metadata::Metadata) = length(metadata.dates)
Base.eltype(metadata::Metadata) = Base.eltype(metadata.dates)

# If only one date, it's a single element array
Base.length(metadata::Metadata{<:AbstractCFDateTime}) = 1

@propagate_inbounds Base.getindex(m::Metadata, i::Int) = Metadata(m.name, m.dates[i],   m.version, m.dir)
@propagate_inbounds Base.first(m::Metadata)            = Metadata(m.name, m.dates[1],   m.version, m.dir)
@propagate_inbounds Base.last(m::Metadata)             = Metadata(m.name, m.dates[end], m.version, m.dir)

@inline function Base.iterate(m::Metadata, i=1)
    if (i % UInt) - 1 < length(m)
        return Metadata(m.name, m.dates[i], m.version, m.dir), i + 1
    else
        return nothing
    end
end

Base.axes(metadata::Metadata{<:AbstractCFDateTime})    = 1
Base.first(metadata::Metadata{<:AbstractCFDateTime})   = metadata
Base.last(metadata::Metadata{<:AbstractCFDateTime})    = metadata
Base.iterate(metadata::Metadata{<:AbstractCFDateTime}) = (metadata, nothing)
Base.iterate(::Metadata{<:AbstractCFDateTime}, ::Any)  = nothing

metadata_path(metadata) = joinpath(metadata.dir, metadata_filename(metadata))


"""
    native_times(metadata; start_time = first(metadata).dates)

Extract the time values from the given metadata and calculates the time difference
from the start time.

Arguments
=========
- `metadata`: The metadata containing the date information.
- `start_time`: The start time for calculating the time difference. Defaults to the first date in the metadata.

Returns
=======
An array of time differences in seconds.
"""
function native_times(metadata; start_time=first(metadata).dates)
    times = zeros(length(metadata))
    for (t, data) in enumerate(metadata)
        date = data.dates
        time = date - start_time
        time = Second(time).value
        times[t] = time
    end

    return times
end

"""
    all_dates(metadata)

Extracts all the dates of the given metadata formatted using the `DateTimeProlepticGregorian` type.
Needs to be extended by any new dataset version.
"""
all_dates(metadata) = all_dates(metadata.version)

# File names of metadata containing multiple dates
metadata_filename(metadata) = [metadata_filename(metadatum) for metadatum in metadata]
