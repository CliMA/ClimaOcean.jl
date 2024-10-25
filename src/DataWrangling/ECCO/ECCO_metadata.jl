using CFTime
using Dates

using Base: @propagate_inbounds

import Oceananigans.Fields: set!
import Base

struct ECCO2Monthly end
struct ECCO2Daily end
struct ECCO4Monthly end

# Metadata holding the ECCO dataset information:
# - `name`: The name of the dataset.
# - `dates`: The dates of the dataset, in a `AbstractCFDateTime` format.
# - `version`: The version of the dataset, could be ECCO2Monthly, ECCO2Daily, or ECCO4Monthly.
# - `path`: The path where the dataset is stored.
struct ECCOMetadata{D, V} 
    name  :: Symbol
    dates :: D
    version :: V
    path :: String
end

Base.show(io::IO, metadata::ECCOMetadata) = 
    print(io, "ECCOMetadata:", '\n',
    "├── field: $(metadata.name)", '\n',
    "├── dates: $(metadata.dates)", '\n',
    "├── version: $(metadata.version)", '\n',
    "└── file path: $(metadata.path)")

"""
    ECCOMetadata(name::Symbol; 
                 date = DateTimeProlepticGregorian(1993, 1, 1), 
                 version = ECCO2Daily(), 
                 path = download_ECCO_cache)

Constructs an `ECCOMetadata` object with the specified parameters.

# Arguments
============
- `name::Symbol`: The name of the metadata.

# Keyword Arguments
===================
- `date`: The date of the metadata (default: DateTimeProlepticGregorian(1993, 1, 1)).
- `version`: The version of the metadata (for the moment the choices are ECCO2Monthly(), ECCO2Daily(), or ECCO4Monthly()).
- `path`: The path to the datafile (default: download_ECCO_cache).
"""
function ECCOMetadata(name::Symbol; 
                      date = DateTimeProlepticGregorian(1993, 1, 1),
                   version = ECCO2Daily(),
                      path = download_ECCO_cache) 
             
    return ECCOMetadata(name, date, version, path) 
end

ECCOMetadata(name::Symbol, date, version=ECCO4Monthly(); path = download_ECCO_cache) = ECCOMetadata(name, date, version, path)

# Treat ECCOMetadata as an array to allow iteration over the dates.
Base.length(metadata::ECCOMetadata)       = length(metadata.dates)
Base.eltype(metadata::ECCOMetadata)       = Base.eltype(metadata.dates)
@propagate_inbounds Base.getindex(m::ECCOMetadata, i::Int) = ECCOMetadata(m.name, m.dates[i],   m.version, m.path)
@propagate_inbounds Base.first(m::ECCOMetadata)            = ECCOMetadata(m.name, m.dates[1],   m.version, m.path)
@propagate_inbounds Base.last(m::ECCOMetadata)             = ECCOMetadata(m.name, m.dates[end], m.version, m.path)

@inline function Base.iterate(m::ECCOMetadata, i=1)
    if (i % UInt) - 1 < length(m)
        return ECCOMetadata(m.name, m.dates[i], m.version, m.path), i + 1
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
    shortname   = short_name(metadata)
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
short_name(data::ECCOMetadata{<:Any, <:ECCO2Daily})   = ECCO2_short_names[data.name]
short_name(data::ECCOMetadata{<:Any, <:ECCO2Monthly}) = ECCO2_short_names[data.name]
short_name(data::ECCOMetadata{<:Any, <:ECCO4Monthly}) = ECCO4_short_names[data.name]

field_location(data::ECCOMetadata) = ECCO_location[data.name]

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
    download_dataset!(metadata::ECCOMetadata;
                      url = urls(metadata))

Download the dataset specified by the given metadata. If the metadata contains a single date, 
the dataset is downloaded directly. If the metadata contains multiple dates, the dataset is 
downloaded for each date individually. 
The data download requires a username and password to be provided in the ECCO_USERNAME and ECCO_PASSWORD
environment variables. This can be done by exporting the environment variables in the shell before running the script,
or by launching julia with 

ECCO_USERNAME=myuser ECCO_PASSWORD=mypasswrd julia 

# Arguments
- `metadata::ECCOMetadata`: The metadata specifying the dataset to be downloaded.
"""
function download_dataset!(metadata::ECCOMetadata;
                           url = urls(metadata))

    username = get(ENV, "ECCO_USERNAME", nothing)
    password = get(ENV, "ECCO_PASSWORD", nothing)
    path     = metadata.path

    for data in metadata
        filename  = metadata_filename(data)
        shortname = short_name(data)

        if !isfile(filename)

            msg = "\n See ClimaOcean.jl/src/ECCO/README.md for instructions."
            isnothing(username) && throw(ArgumentError("Could not find the username for $(url). Please provide a username in the ECCO_USERNAME environment variable." * msg))
            isnothing(password) && throw(ArgumentError("Could not find the username for $(url). Please provide a password in the ECCO_PASSWORD environment variable." * msg))
        
            # Version specific download file url
            if data.version isa ECCO2Monthly || data.version isa ECCO2Daily
                fileurl = joinpath(url, shortname, filename)
            elseif data.version isa ECCO4Monthly
                year    = string(Dates.year(data.dates))
                fileurl = joinpath(url, shortname, year, filename)
            end

            cmd = `wget --http-user=$(username) --http-passwd=$(password)  --directory-prefix=$(path) $(fileurl)`
        
            run(cmd)
        end
    end

    return nothing
end
