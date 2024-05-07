using CFTime
using Dates
import Dates: year, month, day
import Oceananigans.Fields: set!
import Base

struct ECCO2Monthly end
struct ECCO2Daily end
struct ECCO4Monthly end

# Metadata holding the ECCO dataset information:
# - `name`: The name of the dataset.
# - `dates`: The dates of the dataset, in a `AbstractCFDateTime` format.
# - `version`: The version of the dataset, could be ECCO2Monthly, ECCO2Daily, or ECCO4Monthly.
struct ECCOMetadata{D, V} 
    name  :: Symbol
    dates :: D
    version :: V
end

# The default is the ECCO4Monthly dataset at 1992-01-01.
ECCOMetadata(name::Symbol) = ECCOMetadata(name, DateTimeProlepticGregorian(1992, 1, 1), ECCO4Monthly())

# Treat ECCOMetadata as an array to allow iteration over the dates.
Base.getindex(metadata::ECCOMetadata, i::Int) = @inbounds ECCOMetadata(metadata.name, metadata.dates[i], metdata.version)
Base.length(metadata::ECCOMetadata)       = length(metadata.dates)
Base.eltype(metadata::ECCOMetadata)       = Base.eltype(metadata.dates)
Base.first(metadata::ECCOMetadata)        = @inbounds ECCOMetadata(metadata.name, metadata.dates[1], metadata.version)
Base.last(metadata::ECCOMetadata)         = @inbounds ECCOMetadata(metadata.name, metadata.dates[end], metadata.version)
Base.iterate(metadata::ECCOMetadata, i=1) = (@inline; (i % UInt) - 1 < length(metadata) ? (@inbounds ECCOMetadata(metadata.name, metadata.dates[i], metadata.version), i + 1) : nothing)

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
all_ecco_dates(::ECCO4Monthly) = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 1)
all_ecco_dates(::ECCO2Monthly) = DateTimeProlepticGregorian(1992, 1, 1) : Month(1) : DateTimeProlepticGregorian(2023, 12, 1)
all_ecco_dates(::ECCO2Daily)   = DateTimeProlepticGregorian(1992, 1, 4) : Day(1) : DateTimeProlepticGregorian(2023, 12, 31)

# File name generation specific to each Dataset version
function file_name(metadata::ECCOMetadata{<:AbstractCFDateTime, <:ECCO4Monthly})
    shortname   = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return shortname * "_" * yearstr * "_" * monthstr * ".nc"
end

function file_name(metadata::ECCOMetadata{<:AbstractCFDateTime})
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
short_name(data::ECCOMetadata{<:Any, <:ECCO2Daily})   = ecco2_short_names[data.name]
short_name(data::ECCOMetadata{<:Any, <:ECCO2Monthly}) = ecco2_short_names[data.name]
short_name(data::ECCOMetadata{<:Any, <:ECCO4Monthly}) = ecco4_short_names[data.name]

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
urls(::ECCOMetadata{<:Any, <:ECCO2Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/monthly/"
urls(::ECCOMetadata{<:Any, <:ECCO2Daily})   = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/daily/"
urls(::ECCOMetadata{<:Any, <:ECCO4Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/"

"""
    download_dataset!(metadata::ECCOMetadata)

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
                           username = get(ENV, "ECCO_USERNAME", nothing), 
                           password = get(ENV, "ECCO_PASSWORD", nothing),
                           url = urls(metadata))

    isnothing(username) && throw(ArgumentError("Please provide a username in the ECCO_USERNAME environment variable!"))
    isnothing(password) && throw(ArgumentError("Please provide a password in the ECCO_PASSWORD environment variable!"))

    for data in metadata
        filename  = file_name(data)
        shortname = short_name(data)
 
        if !isfile(filename)

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