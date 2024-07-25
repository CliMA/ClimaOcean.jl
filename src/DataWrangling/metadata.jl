
# Metadata holding the ECCO dataset information:
# - `name`: The name of the dataset.
# - `dates`: The dates of the dataset, in a `AbstractCFDateTime` format.
# - `version`: The version of the dataset, could be ECCO2Monthly, ECCO2Daily, or ECCO4Monthly.
struct Metadata{D, V} 
    name  :: Symbol
    dates :: D
    version :: V
end

Base.show(io::IO, metadata::Metadata) = 
    print(io, "Metadata:", '\n',
    "├── field: $(metadata.name)", '\n',
    "├── dates: $(metadata.dates)", '\n',
    "└── data version: $(metadata.version)")

# Treat Metadata as an array to allow iteration over the dates.
Base.getindex(metadata::Metadata, i::Int) = @inbounds Metadata(metadata.name, metadata.dates[i], metadata.version)
Base.length(metadata::Metadata)           = length(metadata.dates)
Base.eltype(metadata::Metadata)           = Base.eltype(metadata.dates)
Base.first(metadata::Metadata)            = @inbounds Metadata(metadata.name, metadata.dates[1], metadata.version)
Base.last(metadata::Metadata)             = @inbounds Metadata(metadata.name, metadata.dates[end], metadata.version)
Base.iterate(metadata::Metadata, i=1)     = (@inline; (i % UInt) - 1 < length(metadata) ? (@inbounds Metadata(metadata.name, metadata.dates[i], metadata.version), i + 1) : nothing)

Base.axes(metadata::Metadata{<:AbstractCFDateTime})    = 1
Base.first(metadata::Metadata{<:AbstractCFDateTime})   = metadata
Base.last(metadata::Metadata{<:AbstractCFDateTime})    = metadata
Base.iterate(metadata::Metadata{<:AbstractCFDateTime}) = (metadata, nothing)
Base.iterate(::Metadata{<:AbstractCFDateTime}, ::Any)  = nothing

Base.size(data::Metadata{<:Any, <:JRA55ThreeHourly}) = (640,  320,  1, length(data.dates))


"""
    native_times(metadata; start_time = metadata.dates[1])

Extracts the time values from the given metadata and calculates the time difference
from the start time.

# Arguments
- `metadata`: The metadata containing the date information.
- `start_time`: The start time for calculating the time difference. Defaults to the first date in the metadata.

# Returns
An array of time differences in seconds.
"""
function native_times(metadata; start_time = first(metadata).dates)
    times = zeros(length(metadata))
    for (t, data) in enumerate(metadata)
        date = data.dates
        time = date - start_time
        time = Second(time).value
        times[t] = time
    end

    return times
end