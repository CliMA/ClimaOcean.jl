module ECCO

export ECCOMetadatum, ECCO_immersed_grid, adjusted_ECCO_tracers, initialize!
export ECCO2Monthly, ECCO4Monthly, ECCO2Daily

using Oceananigans
using ClimaOcean
using NCDatasets
using Dates
using Adapt
using Scratch
using Downloads

using Oceananigans.DistributedComputations: @root

using ClimaOcean.DataWrangling:
    netrc_downloader,
    BoundingBox,
    metadata_path,
    Celsius,
    Metadata,
    Metadatum,
    download_progress

using KernelAbstractions: @kernel, @index

using Dates: year, month, day

import Oceananigans: location

import ClimaOcean.DataWrangling:
    default_download_directory,
    all_dates,
    metadata_filename,
    download_dataset,
    temperature_units,
    dataset_variable_name,
    metaprefix,
    longitude_interfaces,
    latitude_interfaces,
    z_interfaces,
    is_three_dimensional,
    inpainted_metadata_path,
    reversed_vertical_axis,
    default_mask_value,
    reversed_sign

download_ECCO_cache::String = ""
function __init__()
    global download_ECCO_cache = @get_scratch!("ECCO")
end

# Datasets
struct ECCO2Monthly end
struct ECCO2Daily end
struct ECCO4Monthly end
const SomeECCODataset = Union{ECCO2Monthly, ECCO4Monthly, ECCO2Daily}

function default_download_directory(::ECCO2Monthly)
    path = joinpath(download_ECCO_cache, "v2", "monthly")
    return mkpath(path)
end

function default_download_directory(::ECCO2Daily)
    path = joinpath(download_ECCO_cache, "v2", "daily")
    return mkpath(path)
end

function default_download_directory(::ECCO4Monthly)
    path = joinpath(download_ECCO_cache, "v4")
    return mkpath(path)
end

Base.size(::ECCO2Daily, variable)   = (1440, 720, 50)
Base.size(::ECCO2Monthly, variable) = (1440, 720, 50)
Base.size(::ECCO4Monthly, variable) = (720,  360, 50)

temperature_units(::SomeECCODataset) = Celsius()
default_mask_value(::ECCO4Monthly) = 0
reversed_vertical_axis(::SomeECCODataset) = true
reversed_sign(::SomeECCODataset, ::Val{:downwelling_longwave}) = true
reversed_sign(::SomeECCODataset, ::Val{:downwelling_shortwave}) = true

const ECCO2_url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/"
const ECCO4_url = "https://ecco.jpl.nasa.gov/drive/files/Version4/Release4/interp_monthly/"

# The whole range of dates in the different dataset datasets
all_dates(dataset::SomeECCODataset) = all_dates(dataset, nothing)
all_dates(::ECCO4Monthly, variable) = DateTime(1992, 1, 1) : Month(1) : DateTime(2017, 12, 1)
all_dates(::ECCO2Monthly, variable) = DateTime(1992, 1, 1) : Month(1) : DateTime(2024, 12, 1)
all_dates(::ECCO2Daily,   variable) = DateTime(1992, 1, 1) : Day(1)   : DateTime(2024, 12, 31)

longitude_interfaces(::SomeECCODataset) = (-180, 180)
latitude_interfaces(::SomeECCODataset) = (-90, 90)

z_interfaces(::SomeECCODataset) = [
    -6128.75,
    -5683.75,
    -5250.25,
    -4839.75,
    -4452.25,
    -4087.75,
    -3746.25,
    -3427.75,
    -3132.25,
    -2859.75,
    -2610.25,
    -2383.74,
    -2180.13,
    -1999.09,
    -1839.64,
    -1699.66,
    -1575.64,
    -1463.12,
    -1357.68,
    -1255.87,
    -1155.72,
    -1056.53,
    -958.45,
    -862.10,
    -768.43,
    -678.57,
    -593.72,
    -515.09,
    -443.70,
    -380.30,
    -325.30,
    -278.70,
    -240.09,
    -208.72,
    -183.57,
    -163.43,
    -147.11,
    -133.45,
    -121.51,
    -110.59,
    -100.20,
    -90.06,
    -80.01,
    -70.0,
    -60.0,
    -50.0,
    -40.0,
    -30.0,
    -20.0,
    -10.0,
      0.0,
]

ECCO4_dataset_variable_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "EVEL",
    :v_velocity            => "NVEL",
    :free_surface          => "SSH",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_concentration => "SIarea",
    :net_heat_flux         => "oceQnet",
    :sensible_heat_flux    => "EXFhs",
    :latent_heat_flux      => "EXFhl",
    :net_longwave          => "EXFlwnet",
    :downwelling_shortwave => "oceQsw",
    :downwelling_longwave  => "EXFlwdn",    
    :air_temperature       => "EXFatemp",
    :air_specific_humidity => "EXFaqh",
    :sea_level_pressure    => "EXFpress",
    :eastward_wind         => "EXFewind",
    :northward_wind        => "EXFnwind",
    :rain_freshwater_flux  => "EXFpreci",
)

ECCO2_dataset_variable_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "EVEL",
    :v_velocity            => "NVEL",
    :free_surface          => "SSH",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_concentration => "SIarea",
    :net_heat_flux         => "oceQnet",
)

ECCO_location = Dict(
    :temperature           => (Center, Center, Center),
    :salinity              => (Center, Center, Center),
    :u_velocity            => (Face,   Center, Center),
    :v_velocity            => (Center, Face,   Center),
    :free_surface          => (Center, Center, Nothing),
    :sea_ice_thickness     => (Center, Center, Nothing),
    :sea_ice_concentration => (Center, Center, Nothing),
    :net_heat_flux         => (Center, Center, Nothing),
    :sensible_heat_flux    => (Center, Center, Nothing),
    :latent_heat_flux      => (Center, Center, Nothing),
    :net_longwave          => (Center, Center, Nothing),
    :downwelling_longwave  => (Center, Center, Nothing),
    :downwelling_shortwave => (Center, Center, Nothing),
    :air_temperature       => (Center, Center, Nothing),
    :air_specific_humidity => (Center, Center, Nothing),
    :sea_level_pressure    => (Center, Center, Nothing),
    :eastward_wind         => (Center, Center, Nothing),
    :northward_wind        => (Center, Center, Nothing),
    :rain_freshwater_flux  => (Center, Center, Nothing),
)    

const ECCOMetadata{D} = Metadata{<:SomeECCODataset, D}
const ECCOMetadatum   = Metadatum{<:SomeECCODataset}

"""
    ECCOMetadatum(name;
                  date = first_date(ECCO4Monthly(), name),
                  dir = download_ECCO_cache)

An alias to construct a [`Metadatum`](@ref) of `ECCO4Monthly()`.
"""
function ECCOMetadatum(name;
                       date = first_date(ECCO4Monthly(), name),
                       dir = download_ECCO_cache)

    return Metadatum(name; date, dir, dataset=ECCO4Monthly())
end

metaprefix(::ECCOMetadata) = "ECCOMetadata"

# File name generation specific to each dataset
function metadata_filename(metadata::Metadatum{<:ECCO4Monthly})
    shortname = dataset_variable_name(metadata)
    yearstr   = string(Dates.year(metadata.dates))
    monthstr  = string(Dates.month(metadata.dates), pad=2)
    return shortname * "_" * yearstr * "_" * monthstr * ".nc"
end

function metadata_filename(metadata::Metadatum{<:Union{ECCO2Daily, ECCO2Monthly}})
    shortname = dataset_variable_name(metadata)
    yearstr   = string(Dates.year(metadata.dates))
    monthstr  = string(Dates.month(metadata.dates), pad=2)
    postfix   = is_three_dimensional(metadata) ? ".1440x720x50." : ".1440x720."

    if metadata.dataset isa ECCO2Monthly
        return shortname * postfix * yearstr * monthstr * ".nc"
    elseif metadata.dataset isa ECCO2Daily
        daystr = is_three_dimensional(metadata) ? string(Dates.day(metadata.dates), pad=2) : ""
        return shortname * postfix * yearstr * monthstr * daystr * ".nc"
    end
end

# Convenience functions
dataset_variable_name(data::Metadata{<:ECCO2Daily})   = ECCO2_dataset_variable_names[data.name]
dataset_variable_name(data::Metadata{<:ECCO2Monthly}) = ECCO2_dataset_variable_names[data.name]
dataset_variable_name(data::Metadata{<:ECCO4Monthly}) = ECCO4_dataset_variable_names[data.name]
location(data::ECCOMetadata) = ECCO_location[data.name]

is_three_dimensional(data::ECCOMetadata) =
    data.name == :temperature ||
    data.name == :salinity ||
    data.name == :u_velocity ||
    data.name == :v_velocity

# URLs for the ECCO datasets specific to each dataset
metadata_url(m::Metadata{<:ECCO2Monthly}) = ECCO2_url * "monthly/" * dataset_variable_name(m) * "/" * metadata_filename(m)
metadata_url(m::Metadata{<:ECCO2Daily})   = ECCO2_url * "daily/"   * dataset_variable_name(m) * "/" * metadata_filename(m)

function metadata_url(m::Metadata{<:ECCO4Monthly})
    year = string(Dates.year(m.dates))
    return ECCO4_url * dataset_variable_name(m) * "/" * year * "/" * metadata_filename(m)
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

function inpainted_metadata_filename(metadata::ECCOMetadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::ECCOMetadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

include("ECCO_atmosphere.jl")

end # Module
