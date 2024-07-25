using CFTime
using Dates
import Dates: year, month, day
import Oceananigans.Fields: set!
import Base

using ClimaOcean.DataWrangling: Metadata

struct ECCO2Monthly end
struct ECCO2Daily end
struct ECCO4Monthly end

const ECCOMetadata{T} = Union{Metadata{T, <:ECCO4Monthly}, Metadata{T, <:ECCO2Daily}, Metadata{T, <:ECCO2Monthly}} where T

Base.size(data::Metadata{<:Any, <:ECCO2Daily})       = (1440, 720, 50, length(data.dates))
Base.size(data::Metadata{<:Any, <:ECCO2Monthly})     = (1440, 720, 50, length(data.dates))
Base.size(data::Metadata{<:Any, <:ECCO4Monthly})     = (720,  360, 50, length(data.dates))

Base.size(::Metadata{<:AbstractCFDateTime, <:ECCO2Daily})   = (1440, 720, 50, 1)
Base.size(::Metadata{<:AbstractCFDateTime, <:ECCO2Monthly}) = (1440, 720, 50, 1)
Base.size(::Metadata{<:AbstractCFDateTime, <:ECCO4Monthly}) = (720,  360, 50, 1)

# The whole range of dates in the different dataset versions
all_ecco_dates(::ECCO4Monthly) = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 1)
all_ecco_dates(::ECCO2Monthly) = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 1)
all_ecco_dates(::ECCO2Daily)   = DateTimeProlepticGregorian(1992, 1, 4) : Day(1)   : DateTimeProlepticGregorian(2023, 12, 31)

# File name generation specific to each Dataset version
function metadata_filename(metadata::Metadata{<:AbstractCFDateTime, <:ECCO4Monthly})
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
short_name(data::Metadata{<:Any, <:ECCO2Daily})   = ecco2_short_names[data.name]
short_name(data::Metadata{<:Any, <:ECCO2Monthly}) = ecco2_short_names[data.name]
short_name(data::Metadata{<:Any, <:ECCO4Monthly}) = ecco4_short_names[data.name]

field_location(data::ECCOMetadata) = ecco_location[data.name]

variable_is_three_dimensional(data::ECCOMetadata) = 
    data.name == :temperature || 
    data.name == :salinity || 
    data.name == :u_velocity ||
    data.name == :v_velocity

ecco4_short_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "EVEL",
    :v_velocity            => "NVEL",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_area_fraction => "SIarea"
)

ecco2_short_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "UVEL",
    :v_velocity            => "VVEL",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_area_fraction => "SIarea"
)

ecco_location = Dict(
    :temperature           => (Center, Center, Center),
    :salinity              => (Center, Center, Center),
    :sea_ice_thickness     => (Center, Center, Nothing),
    :sea_ice_area_fraction => (Center, Center, Nothing),
    :u_velocity            => (Face,   Center, Center),
    :v_velocity            => (Center, Face,   Center),
)

# URLs for the ECCO datasets specific to each version
urls(::Metadata{<:Any, <:ECCO2Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/monthly/"
urls(::Metadata{<:Any, <:ECCO2Daily})   = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/daily/"
urls(::Metadata{<:Any, <:ECCO4Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/"

"""
    download_dataset!(metadata::Metadata)

Download the dataset specified by the given metadata. If the metadata contains a single date, 
the dataset is downloaded directly. If the metadata contains multiple dates, the dataset is 
downloaded for each date individually. 
The data download requires a username and password to be provided in the ECCO_USERNAME and ECCO_PASSWORD
environment variables. This can be done by exporting the environment variables in the shell before running the script,
or by launching julia with 

ECCO_USERNAME=myuser ECCO_PASSWORD=mypasswrd julia 

# Arguments
- `metadata::Metadata`: The metadata specifying the dataset to be downloaded.
"""
function download_dataset!(metadata::ECCOMetadata;
                           url = urls(metadata))

    username = get(ENV, "ECCO_USERNAME", nothing)
    password = get(ENV, "ECCO_PASSWORD", nothing)

    for data in metadata
        filename  = metadata_filename(data)
        shortname = short_name(data)

        if !isfile(filename)

            isnothing(username) && throw(ArgumentError("Could not find the username for $(url). Please provide a username in the ECCO_USERNAME environment variable."))
            isnothing(password) && throw(ArgumentError("Could not find the username for $(url). Please provide a password in the ECCO_PASSWORD environment variable."))
        
            # Version specific download file url
            if data.version isa ECCO2Monthly || data.version isa ECCO2Daily
                fileurl = joinpath(url, shortname, filename)
            elseif data.version isa ECCO4Monthly
                year    = string(Dates.year(data.dates))
                fileurl = joinpath(url, shortname, year, filename)
            end

            cmd = `wget --http-user=$(username) --http-passwd=$(password) $(fileurl)`
        
            run(cmd)
        end
    end

    return nothing
end