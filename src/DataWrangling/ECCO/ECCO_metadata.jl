using CFTime
using Dates
using Downloads
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: netrc_downloader, metadata_path, AnyDateTime,
                                Celsius, Kelvin
using Oceananigans.DistributedComputations

import Base
import Dates: year, month, day
import Oceananigans.Fields: location
import ClimaOcean.DataWrangling: all_dates, metadata_filename, download_dataset,
                                 default_download_directory, dataset_temperature_units,
                                 short_name, dataset_latitude_extent

struct ECCO2Monthly end
struct ECCO2Daily end
struct ECCO4Monthly end

const ECCOMetadata{D} = Metadata{<:Union{<:ECCO2Monthly, <:ECCO2Daily, <:ECCO4Monthly}, D}
const ECCOMetadatum   = Metadatum{<:Union{<:ECCO2Monthly, <:ECCO2Daily, <:ECCO4Monthly}}

const ECCO2_url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/"
const ECCO4_url = "https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/"

"""
    ECCOMetadatum(name;
                  date = first_date(ECCO4Monthly()),
                  dir = download_ECCO_cache)

An alias to construct a [`Metadatum`](@ref) of `ECCO4Monthly()`.
"""
function ECCOMetadatum(name;
                       date = first_date(ECCO4Monthly()),
                       dir = download_ECCO_cache)

    return Metadatum(name; date, dir, dataset=ECCO4Monthly())
end

default_download_directory(::Union{<:ECCO2Monthly, <:ECCO2Daily, <:ECCO4Monthly}) = download_ECCO_cache

datasetstr(md::ECCOMetadata) = string(md.dataset)

datestr(md::ECCOMetadata) = string(first(md.dates), "--", last(md.dates))
datestr(md::ECCOMetadatum) = string(md.dates)

Base.summary(md::ECCOMetadata) = string("ECCOMetadata{", datasetstr(md), "} of ",
                                        md.name, " for ", datestr(md))

Base.size(data::Metadata{<:ECCO2Daily})   = (1440, 720, 50, length(data.dates))
Base.size(data::Metadata{<:ECCO2Monthly}) = (1440, 720, 50, length(data.dates))
Base.size(data::Metadata{<:ECCO4Monthly}) = (720,  360, 50, length(data.dates))

Base.size(::Metadatum{<:ECCO2Daily})   = (1440, 720, 50, 1)
Base.size(::Metadatum{<:ECCO2Monthly}) = (1440, 720, 50, 1)
Base.size(::Metadatum{<:ECCO4Monthly}) = (720,  360, 50, 1)

# The whole range of dates in the different dataset datasets
all_dates(::ECCO4Monthly, name) = DateTime(1992, 1, 1) : Month(1) : DateTime(2023, 12, 1)
all_dates(::ECCO2Monthly, name) = DateTime(1992, 1, 1) : Month(1) : DateTime(2023, 12, 1)
all_dates(::ECCO2Daily, name)   = DateTime(1992, 1, 4) : Day(1)   : DateTime(2023, 12, 31)

# Fallback, actually, we do not really need the name for ECCO since all
# variables have the same frequency and the same time-range, differently from JRA55
all_dates(dataset::Union{<:ECCO4Monthly, <:ECCO2Monthly, <:ECCO2Daily}) = all_dates(dataset, :temperature)

# File name generation specific to each Dataset dataset
function metadata_filename(metadata::Metadatum{<:ECCO4Monthly})
    shortname = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return shortname * "_" * yearstr * "_" * monthstr * ".nc"
end

function metadata_filename(metadata::Metadatum{<:Union{ECCO2Daily, ECCO2Monthly}})
    shortname   = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    postfix = variable_is_three_dimensional(metadata) ? ".1440x720x50." : ".1440x720."

    if metadata.dataset isa ECCO2Monthly
        return shortname * postfix * yearstr * monthstr * ".nc"
    elseif metadata.dataset isa ECCO2Daily
        daystr = string(Dates.day(metadata.dates), pad=2)
        return shortname * postfix * yearstr * monthstr * daystr * ".nc"
    end
end

# Convenience functions
short_name(data::Metadata{<:ECCO2Daily})   = ECCO2_short_names[data.name]
short_name(data::Metadata{<:ECCO2Monthly}) = ECCO2_short_names[data.name]
short_name(data::Metadata{<:ECCO4Monthly}) = ECCO4_short_names[data.name]

location(data::ECCOMetadata) = ECCO_location[data.name]

dataset_temperature_units(data::ECCOMetadata) = Celsius()

dataset_latitude_extent(data::ECCOMetadatum) = (-90, 90)

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

# URLs for the ECCO datasets specific to each dataset
metadata_url(m::Metadata{<:ECCO2Monthly}) = ECCO2_url * "monthly/" * short_name(m) * "/" * metadata_filename(m)
metadata_url(m::Metadata{<:ECCO2Daily})   = ECCO2_url * "daily/"   * short_name(m) * "/" * metadata_filename(m)

function metadata_url(m::Metadata{<:ECCO4Monthly})
    year = string(Dates.year(m.dates))
    return ECCO4_url * short_name(m) * "/" * year * "/" * metadata_filename(m)
end

function download_dataset(metadata::ECCOMetadata)
    username = get(ENV, "ECCO_USERNAME", nothing)
    password = get(ENV, "ECCO_PASSWORD", nothing)
    dir = metadata.dir

    # Create a temporary directory to store the .netrc file
    # The directory will be deleted after the download is complete
    @root mktempdir(dir) do tmp

        # Write down the username and password in a .netrc file
        downloader = netrc_downloader(username, password, "ecco.jpl.nasa.gov", tmp)
        ntasks = Threads.nthreads()

        asyncmap(metadata; ntasks) do metadatum # Distribute the download among tasks

            fileurl  = metadata_url(metadatum)
            filepath = metadata_path(metadatum)

            if !isfile(filepath)
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
    end

    return nothing
end
