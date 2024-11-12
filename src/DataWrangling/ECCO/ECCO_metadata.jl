using CFTime
using Dates

using Base: @propagate_inbounds

import Oceananigans.Fields: set!, location
import Base

struct ECCO2Monthly end
struct ECCO2Daily end
struct ECCO4Monthly end

"""
    struct ECCOMetadata{D, V}

Metadata information for an ECCO dataset:
- `name`: The name of the dataset.
- `dates`: The dates of the dataset, in an `AbstractCFDateTime` format.
- `version`: The version of the dataset, could be `ECCO2Monthly`, `ECCO2Daily`, or `ECCO4Monthly`.
- `dir`: The directory where the dataset is stored.
"""
struct ECCOMetadata{D, V}
    name  :: Symbol
    dates :: D
    version :: V
    dir :: String
end

Base.show(io::IO, metadata::ECCOMetadata) = 
    print(io, "ECCOMetadata:", '\n',
    "├── name: $(metadata.name)", '\n',
    "├── dates: $(metadata.dates)", '\n',
    "├── version: $(metadata.version)", '\n',
    "└── dir: $(metadata.dir)")

Base.summary(md::ECCOMetadata{<:Any, <:ECCO2Daily})   = "ECCO2Daily $(md.name) metadata ($(first(md.dates))--$(last(md.dates)))"
Base.summary(md::ECCOMetadata{<:Any, <:ECCO2Monthly}) = "ECCO2Monthly $(md.name) metadata ($(first(md.dates))--$(last(md.dates)))"
Base.summary(md::ECCOMetadata{<:Any, <:ECCO4Monthly}) = "ECCO4Monthly $(md.name) metadata ($(first(md.dates))--$(last(md.dates)))"

Base.summary(md::ECCOMetadata{<:AbstractCFDateTime, <:ECCO2Daily})   = "ECCO2Daily $(md.name) metadata at $(md.dates)"
Base.summary(md::ECCOMetadata{<:AbstractCFDateTime, <:ECCO2Monthly}) = "ECCO2Monthly $(md.name) metadata at $(md.dates)"
Base.summary(md::ECCOMetadata{<:AbstractCFDateTime, <:ECCO4Monthly}) = "ECCO4Monthly $(md.name) metadata at $(md.dates)"

"""
    ECCOMetadata(name::Symbol;
                 dates = DateTimeProlepticGregorian(1993, 1, 1),
                 version = ECCO4Monthly(),
                 dir = download_ECCO_cache)

Construct an `ECCOMetadata` object with the specified parameters.

Arguments
=========
- `name::Symbol`: The name of the metadata.

Keyword Arguments
=================
- `dates`: The date(s) of the metadata. Note this can either be a single date,
           representing a snapshot, or a range of dates, representing a time-series.
           Default: `DateTimeProlepticGregorian(1993, 1, 1)`.

- `version`: The data version. Supported versions are `ECCO2Monthly()`, `ECCO2Daily()`,
             or `ECCO4Monthly()`.

- `dir`: The directory of the data file. Default: `download_ECCO_cache`.
"""
function ECCOMetadata(name::Symbol; 
                      dates = DateTimeProlepticGregorian(1993, 1, 1),
                      version = ECCO4Monthly(),
                      dir = download_ECCO_cache)

    return ECCOMetadata(name, dates, version, dir)
end

# Treat ECCOMetadata as an array to allow iteration over the dates.
Base.length(metadata::ECCOMetadata) = length(metadata.dates)
Base.eltype(metadata::ECCOMetadata) = Base.eltype(metadata.dates)

@propagate_inbounds Base.getindex(m::ECCOMetadata, i::Int) = ECCOMetadata(m.name, m.dates[i],   m.version, m.dir)
@propagate_inbounds Base.first(m::ECCOMetadata)            = ECCOMetadata(m.name, m.dates[1],   m.version, m.dir)
@propagate_inbounds Base.last(m::ECCOMetadata)             = ECCOMetadata(m.name, m.dates[end], m.version, m.dir)

@inline function Base.iterate(m::ECCOMetadata, i=1)
    if (i % UInt) - 1 < length(m)
        return ECCOMetadata(m.name, m.dates[i], m.version, m.dir), i + 1
    else
        return nothing
    end
end

Base.axes(metadata::ECCOMetadata{<:AbstractCFDateTime})    = 1
Base.first(metadata::ECCOMetadata{<:AbstractCFDateTime})   = metadata
Base.last(metadata::ECCOMetadata{<:AbstractCFDateTime})    = metadata
Base.iterate(metadata::ECCOMetadata{<:AbstractCFDateTime}) = (metadata, nothing)
Base.iterate(::ECCOMetadata{<:AbstractCFDateTime}, ::Any)  = nothing

Base.size(data::ECCOMetadata{<:Any, <:ECCO2Daily})   = (1440, 720, 50, length(data.dates))
Base.size(data::ECCOMetadata{<:Any, <:ECCO2Monthly}) = (1440, 720, 50, length(data.dates))
Base.size(data::ECCOMetadata{<:Any, <:ECCO4Monthly}) = (720,  360, 50, length(data.dates))

Base.size(::ECCOMetadata{<:AbstractCFDateTime, <:ECCO2Daily})   = (1440, 720, 50, 1)
Base.size(::ECCOMetadata{<:AbstractCFDateTime, <:ECCO2Monthly}) = (1440, 720, 50, 1)
Base.size(::ECCOMetadata{<:AbstractCFDateTime, <:ECCO4Monthly}) = (720,  360, 50, 1)

# The whole range of dates in the different dataset versions
all_ECCO_dates(::ECCO4Monthly) = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 1)
all_ECCO_dates(::ECCO2Monthly) = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 1)
all_ECCO_dates(::ECCO2Daily)   = DateTimeProlepticGregorian(1992, 1, 4) : Day(1)   : DateTimeProlepticGregorian(2023, 12, 31)

# File name generation specific to each Dataset version
function metadata_filename(metadata::ECCOMetadata{<:AbstractCFDateTime, <:ECCO4Monthly})
    shortname = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return shortname * "_" * yearstr * "_" * monthstr * ".nc"
end

function metadata_filename(metadata::ECCOMetadata{<:AbstractCFDateTime})
    shortname   = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    postfix = variable_is_three_dimensional(metadata) ? ".1440x720x50." : ".1440x720."

    if metadata.version isa ECCO2Monthly 
        return shortname * postfix * yearstr * monthstr * ".nc"
    elseif metadata.version isa ECCO2Daily
        daystr = string(Dates.day(metadata.dates), pad=2)
        return shortname * postfix * yearstr * monthstr * daystr * ".nc"
    end
end

# Convenience functions
metadata_path(metadata) = joinpath(metadata.dir, metadata_filename(metadata))
short_name(data::ECCOMetadata{<:Any, <:ECCO2Daily})   = ECCO2_short_names[data.name]
short_name(data::ECCOMetadata{<:Any, <:ECCO2Monthly}) = ECCO2_short_names[data.name]
short_name(data::ECCOMetadata{<:Any, <:ECCO4Monthly}) = ECCO4_short_names[data.name]

metadata_url(prefix, m::ECCOMetadata{<:Any, <:ECCO2Daily}) = joinpath(prefix, short_name(m), metadata_filename(m))
metadata_url(prefix, m::ECCOMetadata{<:Any, <:ECCO2Monthly}) = joinpath(prefix, short_name(m), metadata_filename(m))

function metadata_url(prefix, m::ECCOMetadata{<:Any, <:ECCO4Monthly})
    year = string(Dates.year(m.dates))
    return joinpath(prefix, short_name(m), year, metadata_filename(m))
end

location(data::ECCOMetadata) = ECCO_location[data.name]

variable_is_three_dimensional(data::ECCOMetadata) =
    data.name == :temperature || 
    data.name == :salinity || 
    data.name == :u_velocity ||
    data.name == :v_velocity

ECCO4_short_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "EVEL",
    :v_velocity            => "NVEL",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_area_fraction => "SIarea"
)

ECCO2_short_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "UVEL",
    :v_velocity            => "VVEL",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_area_fraction => "SIarea"
)

ECCO_location = Dict(
    :temperature           => (Center, Center, Center),
    :salinity              => (Center, Center, Center),
    :sea_ice_thickness     => (Center, Center, Nothing),
    :sea_ice_area_fraction => (Center, Center, Nothing),
    :u_velocity            => (Face,   Center, Center),
    :v_velocity            => (Center, Face,   Center),
)

# URLs for the ECCO datasets specific to each version
urls(::ECCOMetadata{<:Any, <:ECCO2Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/monthly/"
urls(::ECCOMetadata{<:Any, <:ECCO2Daily})   = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/daily/"
urls(::ECCOMetadata{<:Any, <:ECCO4Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/"

"""
    download_dataset!(metadata::ECCOMetadata; url = urls(metadata))

Download the dataset specified by `ECCOMetadata`. If `ECCOMetadata.dates` is a single date, 
the dataset is downloaded directly. If `ECCOMetadata.dates` is a vector of dates, each date
is downloaded individually.
The data download requires a username and password to be provided in the `ECCO_USERNAME` and
`ECCO_PASSWORD` environment variables. This can be done by exporting the environment variables
in the shell before running the script, or by launching julia with 

```
ECCO_USERNAME=myusername ECCO_PASSWORD=mypassword julia 
```

Arguments
=========
- `metadata::ECCOMetadata`: The metadata specifying the dataset to be downloaded.
"""
function download_dataset!(metadata::ECCOMetadata; url = urls(metadata))
    username = get(ENV, "ECCO_USERNAME", nothing)
    password = get(ENV, "ECCO_PASSWORD", nothing)
    dir = metadata.dir

    for metadatum in metadata
        filename = metadata_filename(metadatum)
        filepath = metadata_path(metadatum)

        if !isfile(filepath)
            instructions_msg = "\n See ClimaOcean.jl/src/ECCO/README.md for instructions."
            if isnothing(username)
                msg = "Could not find the ECCO_PASSWORD environment variable. \
                       See ClimaOcean.jl/src/ECCO/README.md for instructions on obtaining \
                       and setting your ECCO_USERNAME and ECCO_PASSWORD."
                throw(ArgumentError(msg))
            elseif isnothing(password)
                msg = "Could not find the ECCO_PASSWORD environment variable. \
                       See ClimaOcean.jl/src/ECCO/README.md for instructions on obtaining \
                       and setting your ECCO_USERNAME and ECCO_PASSWORD."
                throw(ArgumentError(msg))
            end

            fileurl = metadata_url(url, metadatum) 
            cmd = `wget --http-user=$(username) --http-passwd=$(password) --directory-prefix=$dir $fileurl`
            run(cmd)
        end
    end

    return nothing
end
