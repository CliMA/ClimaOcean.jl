using CFTime
using Dates
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: netrc_downloader, metadata_path, AnyDateTime

import Dates: year, month, day
using Downloads

import Oceananigans.Fields: set!, location
import Base
import ClimaOcean.DataWrangling: all_dates, metadata_filename, default_download_folder

struct ECCO2Monthly end
struct ECCO2Daily end
struct ECCO4Monthly end

const ECCOMetadata{D} = Metadata{D, <:Union{<:ECCO2Monthly, <:ECCO2Daily, <:ECCO4Monthly}} where {D}
const ECCOMetadatum   = ECCOMetadata{<:AnyDateTime}

default_download_folder(::Union{<:ECCO2Monthly, <:ECCO2Daily, <:ECCO4Monthly}) = download_ECCO_cache

versionstr(md::ECCOMetadata) = string(md.version)

datestr(md::ECCOMetadata) = string(first(md.dates), "--", last(md.dates))
datestr(md::ECCOMetadatum) = string(md.dates)

Base.summary(md::ECCOMetadata) = string("ECCOMetadata{", versionstr(md), "} of ",
                                        md.name, " for ", datestr(md))

Base.size(data::Metadata{<:Any, <:ECCO2Daily})   = (1440, 720, 50, length(data.dates))
Base.size(data::Metadata{<:Any, <:ECCO2Monthly}) = (1440, 720, 50, length(data.dates))
Base.size(data::Metadata{<:Any, <:ECCO4Monthly}) = (720,  360, 50, length(data.dates))

Base.size(::Metadata{<:AnyDateTime, <:ECCO2Daily})   = (1440, 720, 50, 1)
Base.size(::Metadata{<:AnyDateTime, <:ECCO2Monthly}) = (1440, 720, 50, 1)
Base.size(::Metadata{<:AnyDateTime, <:ECCO4Monthly}) = (720,  360, 50, 1)

# The whole range of dates in the different dataset versions
all_dates(::ECCO4Monthly, name) = DateTime(1992, 1, 1) : Month(1) : DateTime(2023, 12, 1)
all_dates(::ECCO2Monthly, name) = DateTime(1992, 1, 1) : Month(1) : DateTime(2023, 12, 1)
all_dates(::ECCO2Daily, name)   = DateTime(1992, 1, 4) : Day(1)   : DateTime(2023, 12, 31)

# Fallback, actually, we do not really need the name for ECCO since all 
# variables have the same frequency and the same time-range, differently from JRA55
all_dates(version::Union{<:ECCO4Monthly, <:ECCO2Monthly, <:ECCO2Daily}) = all_dates(version, :temperature)

# File name generation specific to each Dataset version
function metadata_filename(metadata::Metadata{<:AnyDateTime, <:ECCO4Monthly})
    shortname = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return shortname * "_" * yearstr * "_" * monthstr * ".nc"
end

function metadata_filename(metadata::Metadata{<:AnyDateTime, <:Union{ECCO2Daily, ECCO2Monthly}})
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
short_name(data::Metadata{<:Any, <:ECCO2Daily})   = ECCO2_short_names[data.name]
short_name(data::Metadata{<:Any, <:ECCO2Monthly}) = ECCO2_short_names[data.name]
short_name(data::Metadata{<:Any, <:ECCO4Monthly}) = ECCO4_short_names[data.name]

metadata_url(prefix, m::Metadata{<:Any, <:ECCO2Daily})   = prefix * "/" * short_name(m) * "/" * metadata_filename(m)
metadata_url(prefix, m::Metadata{<:Any, <:ECCO2Monthly}) = prefix * "/" * short_name(m) * "/" * metadata_filename(m)

function metadata_url(prefix, m::Metadata{<:Any, <:ECCO4Monthly})
    year = string(Dates.year(m.dates))
    return prefix * "/" * short_name(m) * "/" * year * "/" * metadata_filename(m)
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
    :free_surface          => "SSH",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_concentration => "SIarea",
    :net_heat_flux         => "oceQnet"
)

ECCO2_short_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "UVEL",
    :v_velocity            => "VVEL",
    :free_surface          => "SSH",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_concentration => "SIarea",
    :net_heat_flux         => "oceQnet"
)

ECCO_location = Dict(
    :temperature           => (Center, Center, Center),
    :salinity              => (Center, Center, Center),
    :free_surface          => (Center, Center, Nothing),
    :sea_ice_thickness     => (Center, Center, Nothing),
    :sea_ice_concentration => (Center, Center, Nothing),
    :net_heat_flux         => (Center, Center, Nothing),
    :u_velocity            => (Face,   Center, Center),
    :v_velocity            => (Center, Face,   Center),
)

# URLs for the ECCO datasets specific to each version
urls(::Metadata{<:Any, <:ECCO2Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/monthly"
urls(::Metadata{<:Any, <:ECCO2Daily})   = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/daily"
urls(::Metadata{<:Any, <:ECCO4Monthly}) = "https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly"

"""
    download_dataset(metadata::ECCOMetadata; url = urls(metadata))

Download the dataset specified by `ECCOMetadata`. If `ECCOMetadata.dates` is a single date,
the dataset is downloaded directly. If `ECCOMetadata.dates` is a vector of dates, each date
is downloaded individually.

The data download requires a username and password to be provided in the `ECCO_USERNAME` and
`ECCO_PASSWORD` environment variables. This can be done by exporting the environment variables
in the shell before running the script, or by launching julia with

```
ECCO_USERNAME=myusername ECCO_PASSWORD=mypassword julia
```

or by invoking

```julia
julia> ENV["ECCO_USERNAME"] = "myusername"

julia> ENV["ECCO_PASSWORD"] = "mypassword"
```

within julia.


Arguments
=========
- `metadata::ECCOMetadata`: The metadata specifying the dataset to be downloaded.
"""
function download_dataset(metadata::ECCOMetadata; url = urls(metadata))
    username = get(ENV, "ECCO_USERNAME", nothing)
    password = get(ENV, "ECCO_PASSWORD", nothing)
    dir = metadata.dir

    # Create a temporary directory to store the .netrc file
    # The directory will be deleted after the download is complete
    @root mktempdir(dir) do tmp

        # Write down the username and password in a .netrc file
        downloader = netrc_downloader(username, password, "ecco.jpl.nasa.gov", tmp)
        ntasks = Threads.nthreads()
       
        missing_files = false

        asyncmap(metadata; ntasks) do metadatum # Distribute the download among tasks

            fileurl  = metadata_url(url, metadatum) 
            filepath = metadata_path(metadatum)

            if !isfile(filepath)
                missing_files = true
                instructions_msg = "\n See ClimaOcean.jl/src/DataWrangling/ECCO/README.md for instructions."
                if isnothing(username)
                    msg = "Could not find the ECCO_PASSWORD environment variable. \
                            See ClimaOcean.jl/src/DataWrangling/ECCO/README.md for instructions on obtaining \
                            and setting your ECCO_USERNAME and ECCO_PASSWORD." * instructions_msg
                    throw(ArgumentError(msg))
                elseif isnothing(password)
                    msg = "Could not find the ECCO_PASSWORD environment variable. \
                            See ClimaOcean.jl/src/DataWrangling/ECCO/README.md for instructions on obtaining \
                            and setting your ECCO_USERNAME and ECCO_PASSWORD." * instructions_msg
                    throw(ArgumentError(msg))
                end
                @info "Downloading ECCO data: $(metadatum.name) in $(metadatum.dir)..."
                Downloads.download(fileurl, filepath; downloader, progress=download_progress)
            end
        end
        
	    if !missing_files
            @info "Note: ECCO $(metadata.name) data is in $(metadata.dir)."
	    end
    end

    return nothing
end
