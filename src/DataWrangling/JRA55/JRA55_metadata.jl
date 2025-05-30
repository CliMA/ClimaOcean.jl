using CFTime
using Dates
using Downloads

using Oceananigans.DistributedComputations

using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: Metadata, metadata_path, download_progress, AnyDateTime

import Dates: year, month, day
import Oceananigans.Fields: set!
import Base

import Oceananigans.Fields: set!, location
import ClimaOcean.DataWrangling: all_dates, metadata_filename, download_dataset, default_download_directory, available_variables

struct MultiYearJRA55 end
struct RepeatYearJRA55 end

const JRA55Metadata{D} = Metadata{<:Union{<:MultiYearJRA55, <:RepeatYearJRA55}, D}
const JRA55Metadatum   = Metadatum{<:Union{<:MultiYearJRA55, <:RepeatYearJRA55}}

default_download_directory(::Union{<:MultiYearJRA55, <:RepeatYearJRA55}) = download_JRA55_cache

Base.size(data::JRA55Metadata) = (640, 320, length(data.dates))
Base.size(::JRA55Metadatum)    = (640, 320, 1)

# JRA55 is a spatially 2D dataset
is_three_dimensional(data::JRA55Metadata) = false

# The whole range of dates in the different dataset datasets
# NOTE! rivers and icebergs have a different frequency! (typical JRA55 data is three-hourly while rivers and icebergs are daily)
function all_dates(::RepeatYearJRA55, name)
    if name == :river_freshwater_flux || name == :iceberg_freshwater_flux
        return DateTime(1990, 1, 1) : Day(1) : DateTime(1990, 12, 31)
    else
        return DateTime(1990, 1, 1) : Hour(3) : DateTime(1990, 12, 31, 23, 59, 59)
    end
end

all_dates(::MultiYearJRA55, name) = JRA55_multiple_year_dates[name]

# Fallback, if we not provide the name, take the highest frequency
all_dates(dataset::Union{<:MultiYearJRA55, <:RepeatYearJRA55}) = all_dates(dataset, :temperature)

# Valid for all JRA55 datasets
function JRA55_time_indices(dataset, dates, name)
    all_JRA55_dates = all_dates(dataset, name)
    indices = Int[]

    for date in dates
        index = findfirst(x -> x == date, all_JRA55_dates)
        !isnothing(index) && push!(indices, index)
    end

    return indices
end

# File name generation specific to each Dataset dataset
# Note that `RepeatYearJRA55` has only one file associated, so we can define
# the filename directly for the whole `Metadata` object, independent of the `dates`
function metadata_filename(metadata::Metadata{<:RepeatYearJRA55}) # No difference
    shortname = dataset_variable_name(metadata)
    return "RYF." * shortname * ".1990_1991.nc"
end

function metadata_filename(metadata::Metadatum{<:MultiYearJRA55})
    # fix the filename
    shortname = dataset_variable_name(metadata)
    year      = Dates.year(metadata.dates)
    suffix    = "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_"

    end_date = last(JRA55_multiple_year_dates[metadata.name])
    end_hour = Hour(end_date)

    if end_hour == Hour(0)
        dates = "$(year)0101-$(year)1231"
    elseif end_hour == Hour(22)
        dates = "$(year)01010130-$(year)12312230"
    else
        dates = "$(year)01010000-$(year)12312100"
    end

    return shortname * suffix * dates * ".nc"
end

# Convenience functions
dataset_variable_name(data::JRA55Metadata) = JRA55_dataset_variable_names[data.name]
location(::JRA55Metadata) = (Center, Center, Center)

available_variables(::MultiYearJRA55)  = JRA55_variable_names
available_variables(::RepeatYearJRA55) = JRA55_variable_names

# A list of all variables provided in the JRA55 dataset:
JRA55_variable_names = (:river_freshwater_flux,
                        :rain_freshwater_flux,
                        :snow_freshwater_flux,
                        :iceberg_freshwater_flux,
                        :specific_humidity,
                        :sea_level_pressure,
                        :downwelling_longwave_radiation,
                        :downwelling_shortwave_radiation,
                        :temperature,
                        :eastward_velocity,
                        :northward_velocity)

JRA55_dataset_variable_names = Dict(
    :river_freshwater_flux           => "friver",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "prra",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "prsn",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "licalvf",  # Freshwater flux from calving icebergs
    :specific_humidity               => "huss",     # Surface specific humidity
    :sea_level_pressure              => "psl",      # Sea level pressure
    :downwelling_longwave_radiation  => "rlds",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "rsds",     # Downwelling shortwave radiation
    :temperature                     => "tas",      # Near-surface air temperature
    :eastward_velocity               => "uas",      # Eastward near-surface wind
    :northward_velocity              => "vas",      # Northward near-surface wind
)

JRA55_multiple_year_url = "https://esgf-data2.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/OMIP/MRI/MRI-JRA55-do-1-5-0/"

JRA55_multiple_year_prefix = Dict(
    :river_freshwater_flux           => "land/day",
    :rain_freshwater_flux            => "atmos/3hr",
    :snow_freshwater_flux            => "atmos/3hr",
    :iceberg_freshwater_flux         => "landIce/day",
    :specific_humidity               => "atmos/3hrPt",
    :sea_level_pressure              => "atmos/3hrPt",
    :downwelling_longwave_radiation  => "atmos/3hr",
    :downwelling_shortwave_radiation => "atmos/3hr",
    :temperature                     => "atmos/3hrPt",
    :eastward_velocity               => "atmos/3hrPt",
    :northward_velocity              => "atmos/3hrPt",
)

JRA55_multiple_year_dates = Dict(
    :river_freshwater_flux           => DateTime(1958, 1, 1)        : Day(1)  : DateTime(2019, 12, 31),
    :rain_freshwater_flux            => DateTime(1958, 1, 1, 1, 30) : Hour(3) : DateTime(2019, 12, 31, 22, 30),
    :snow_freshwater_flux            => DateTime(1958, 1, 1, 1, 30) : Hour(3) : DateTime(2019, 12, 31, 22, 30),
    :iceberg_freshwater_flux         => DateTime(1958, 1, 1)        : Day(1)  : DateTime(2019, 12, 31),
    :specific_humidity               => DateTime(1958, 1, 1)        : Hour(3) : DateTime(2019, 12, 31, 21),
    :sea_level_pressure              => DateTime(1958, 1, 1)        : Hour(3) : DateTime(2019, 12, 31, 21),
    :downwelling_longwave_radiation  => DateTime(1958, 1, 1, 1, 30) : Hour(3) : DateTime(2019, 12, 31, 22, 30),
    :downwelling_shortwave_radiation => DateTime(1958, 1, 1, 1, 30) : Hour(3) : DateTime(2019, 12, 31, 22, 30),
    :temperature                     => DateTime(1958, 1, 1)        : Hour(3) : DateTime(2019, 12, 31, 21),
    :eastward_velocity               => DateTime(1958, 1, 1)        : Hour(3) : DateTime(2019, 12, 31, 21),
    :northward_velocity              => DateTime(1958, 1, 1)        : Hour(3) : DateTime(2019, 12, 31, 21)
)

JRA55_repeat_year_urls = Dict(
    :shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                            "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :river_freshwater_flux => "https://www.dropbox.com/scl/fi/21ggl4p74k4zvbf04nb67/" *
                              "RYF.friver.1990_1991.nc?rlkey=ny2qcjkk1cfijmwyqxsfm68fz&dl=0",

    :rain_freshwater_flux => "https://www.dropbox.com/scl/fi/5icl1gbd7f5hvyn656kjq/" *
                             "RYF.prra.1990_1991.nc?rlkey=iifyjm4ppwyd8ztcek4dtx0k8&dl=0",

    :snow_freshwater_flux => "https://www.dropbox.com/scl/fi/1r4ajjzb3643z93ads4x4/" *
                             "RYF.prsn.1990_1991.nc?rlkey=auyqpwn060cvy4w01a2yskfah&dl=0",

    :iceberg_freshwater_flux => "https://www.dropbox.com/scl/fi/44nc5y27ohvif7lkvpyv0/" *
                                "RYF.licalvf.1990_1991.nc?rlkey=w7rqu48y2baw1efmgrnmym0jk&dl=0",

    :specific_humidity => "https://www.dropbox.com/scl/fi/66z6ymfr4ghkynizydc29/" *
                          "RYF.huss.1990_1991.nc?rlkey=107yq04aew8lrmfyorj68v4td&dl=0",

    :sea_level_pressure => "https://www.dropbox.com/scl/fi/0fk332027oru1iiseykgp/" *
                           "RYF.psl.1990_1991.nc?rlkey=4xpr9uah741483aukok6d7ctt&dl=0",

    :downwelling_longwave_radiation  => "https://www.dropbox.com/scl/fi/y6r62szkirrivua5nqq61/" *
                                        "RYF.rlds.1990_1991.nc?rlkey=wt9yq3cyrvs2rbowoirf4nkum&dl=0",

    :downwelling_shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                                        "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :temperature => "https://www.dropbox.com/scl/fi/fpl0npwi476w635g6lke9/" *
                    "RYF.tas.1990_1991.nc?rlkey=0skb9pe6lgbfbiaoybe7m945s&dl=0",

    :eastward_velocity => "https://www.dropbox.com/scl/fi/86wetpqla2x97isp8092g/" *
                          "RYF.uas.1990_1991.nc?rlkey=rcaf18sh1yz0v9g4hjm1249j0&dl=0",

    :northward_velocity => "https://www.dropbox.com/scl/fi/d38sflo9ddljstd5jwgml/" *
                           "RYF.vas.1990_1991.nc?rlkey=f9y3e57kx8xrb40gbstarf0x6&dl=0",
)

metadata_url(metadata::Metadata{<:RepeatYearJRA55}) = JRA55_repeat_year_urls[metadata.name]

function metadata_url(m::Metadata{<:MultiYearJRA55})
    prefix = JRA55_multiple_year_prefix[m.name]
    return JRA55_multiple_year_url * prefix * "/" * dataset_variable_name(m) * "/gr/v20200916/" * metadata_filename(m)
end

function download_dataset(metadata::JRA55Metadata)

    @root for metadatum in metadata

        fileurl  = metadata_url(metadatum)
        filepath = metadata_path(metadatum)

        if !isfile(filepath)
            Downloads.download(fileurl, filepath; progress=download_progress)
        end
    end

    return nothing
end
