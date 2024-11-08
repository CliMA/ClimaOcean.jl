

"""
    Metadata{D, V} 

Metadata holding the dataset information:
- `name`: The name of the dataset.
- `dates`: The dates of the dataset, in a `AbstractCFDateTime` format.
- `version`: The version of the dataset, could be ECCO2Monthly, ECCO2Daily, or ECCO4Monthly.
- `dir`: The directory where the dataset is stored.
"""
struct Metadata{D, V} 
    name  :: Symbol
    dates :: D
    version :: V
    dir :: String
end

Base.show(io::IO, metadata::Metadata) = 
    print(io, "Metadata:", '\n',
    "├── field: $(metadata.name)", '\n',
    "├── dates: $(metadata.dates)", '\n',
    "└── data version: $(metadata.version)")


# Treat Metadata as an array to allow iteration over the dates.
Base.length(metadata::Metadata) = length(metadata.dates)
Base.eltype(metadata::Metadata) = Base.eltype(metadata.dates)
@propagate_inbounds Base.getindex(m::Metadata, i::Int) = Metadata(m.name, m.dates[i],   m.version, m.dir)
@propagate_inbounds Base.first(m::Metadata)            = Metadata(m.name, m.dates[1],   m.version, m.dir)
@propagate_inbounds Base.last(m::Metadata)             = Metadata(m.name, m.dates[end], m.version, m.dir)

@inline function Base.iterate(m::Metadata, i=1)
    if (i % UInt) - 1 < length(m)
        return ECCOMetadata(m.name, m.dates[i], m.version, m.dir), i + 1
    else
        return nothing
    end
end

Base.axes(metadata::Metadata{<:AbstractCFDateTime})    = 1
Base.first(metadata::Metadata{<:AbstractCFDateTime})   = metadata
Base.last(metadata::Metadata{<:AbstractCFDateTime})    = metadata
Base.iterate(metadata::Metadata{<:AbstractCFDateTime}) = (metadata, nothing)
Base.iterate(::Metadata{<:AbstractCFDateTime}, ::Any)  = nothing

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