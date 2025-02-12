using CFTime
using Dates
using Downloads

using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: Metadata, metadata_path, download_progress

import Dates: year, month, day
import Oceananigans.Fields: set!
import Base
import ClimaOcean.DataWrangling: all_dates, metadata_filename

import Oceananigans.Fields: set!, location
import ClimaOcean.DataWrangling: all_dates, metadata_filename

struct JRA55MultipleYears end
struct JRA55RepeatYear end

const JRA55Metadata{T, V} = Metadata{T, V} where {T, V<:Union{<:JRA55MultipleYears, <:JRA55RepeatYear}}

Base.size(data::JRA55Metadata) = (640, 320, length(data.dates))
Base.size(::JRA55Metadata{<:AbstractCFDateTime}) = (640, 320, 1)

# The whole range of dates in the different dataset versions
all_dates(::JRA55RepeatYear)    = DateTimeProlepticGregorian(1990, 1, 1) : Hour(3) : DateTimeProlepticGregorian(1991, 1, 1)
all_dates(::JRA55MultipleYears) = DateTimeProlepticGregorian(1958, 1, 1) : Hour(3) : DateTimeProlepticGregorian(2021, 1, 1)

function JRA55_time_indices(version, dates)
    all_JRA55_dates = all_dates(version)
    indices = Int[]
    
    for date in dates
        index = findfirst(x -> x == date, all_JRA55_dates)
        push!(indices, index)
    end

    return indices
end

# File name generation specific to each Dataset version
function metadata_filename(metadata::Metadata{<:Any, <:JRA55RepeatYear}) # No difference 
    shortname = short_name(metadata)
    return "RYF." * shortname * ".1990_1991.nc"
end

function metadata_filename(metadata::Metadata{<:AbstractCFDateTime, <:JRA55MultipleYears})
    # fix the filename
    shortname = short_name(metadata)
    year      = Dates.year(metadata.dates)
    suffix    = "_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_"
    dates     = "($year)01010130-($year)12312330.nc"
    return shortname * suffix * dates * ".nc"
end

# Convenience functions
short_name(data::JRA55Metadata) = JRA55_short_names[data.name]
location(::JRA55Metadata)   = (Center, Center, Center)

# A list of all variables provided in the JRA55 dataset:
JRA55_variable_names = (:river_freshwater_flux,
                        :rain_freshwater_flux,
                        :snow_freshwater_flux,
                        :iceberg_freshwater_flux,
                        :specific_humidity,
                        :sea_level_pressure,
                        :relative_humidity,
                        :downwelling_longwave_radiation,
                        :downwelling_shortwave_radiation,
                        :temperature,
                        :eastward_velocity,
                        :northward_velocity)

JRA55_short_names = Dict(
    :river_freshwater_flux           => "friver",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "prra",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "prsn",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "licalvf",  # Freshwater flux from calving icebergs
    :specific_humidity               => "huss",     # Surface specific humidity
    :sea_level_pressure              => "psl",      # Sea level pressure
    :relative_humidity               => "rhuss",    # Surface relative humidity
    :downwelling_longwave_radiation  => "rlds",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "rsds",     # Downwelling shortwave radiation
    :temperature                     => "tas",      # Near-surface air temperature
    :eastward_velocity               => "uas",      # Eastward near-surface wind
    :northward_velocity              => "vas",      # Northward near-surface wind
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

    :relative_humidity => "https://www.dropbox.com/scl/fi/1agwsp0lzvntuyf8bm9la/" *
                          "RYF.rhuss.1990_1991.nc?rlkey=8cd0vs7iy1rw58b9pc9t68gtz&dl=0",

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

variable_is_three_dimensional(data::JRA55Metadata) = false

urls(metadata::Metadata{<:Any, <:JRA55RepeatYear}) = JRA55_repeat_year_urls[metadata.name]  
function urls(metadata::Metadata{<:Any, <:JRA55MultipleYears})
    shotname = short_name(metadata)
    return "https://esgf-data2.llnl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/OMIP/MRI/MRI-JRA55-do-1-5-0/atmos/$(shortname)/prra/gr/v20200916/"
end

metadata_url(prefix, m::Metadata{<:Any, <:JRA55RepeatYear}) = prefix # No specific name for this url

# TODO: This will need to change when we add a method for JRA55MultipleYears
function download_dataset!(metadata::JRA55Metadata; url = urls(metadata))

    asyncmap(metadata, ntasks=10) do metadatum # Distribute the download among tasks

        fileurl  = metadata_url(url, metadatum) 
        filepath = metadata_path(metadatum)

        if !isfile(filepath)
            Downloads.download(fileurl, filepath; progress=download_progress)
        end
    end

    return nothing
end