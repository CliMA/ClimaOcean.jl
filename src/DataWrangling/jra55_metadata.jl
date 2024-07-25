using CFTime
using Dates
import Dates: year, month, day
import Oceananigans.Fields: set!
import Base

using ClimaOcean.DataWrangling: Metadata

struct JRA55MultipleYears end
struct JRA55RepeatYear end

const JRA55Metadata{T} = Union{Metadata{T, <:JRA55MultipleYears}, Metadata{<:JRA55RepeatYear}} where T

Base.size(data::JRA55Metadata) = (640, 320, length(data.dates))
Base.size(::JRA55Metadata{<:AbstractCFDateTime}) = (640, 320, 1)

# File name generation specific to each Dataset version
function metadata_filename(metadata::Metadata{<:AbstractCFDateTime, <:JRA55MultipleYears})
    # fix the filename

    return filename
end

# File name generation specific to each Dataset version
function metadata_filename(metadata::Metadata{<:AbstractCFDateTime, <:JRA55RepeatYear})
    # fix the filename

    return filename
end

# Convenience functions
short_name(data::JRA55Metadata)     = jra55_short_names[data.name]
field_location(data::JRA55Metadata) = jra55_location[data.name]

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

filenames = Dict(
    :river_freshwater_flux           => "RYF.friver.1990_1991.nc",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "RYF.prra.1990_1991.nc",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "RYF.prsn.1990_1991.nc",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "RYF.licalvf.1990_1991.nc",  # Freshwater flux from calving icebergs
    :specific_humidity               => "RYF.huss.1990_1991.nc",     # Surface specific humidity
    :sea_level_pressure              => "RYF.psl.1990_1991.nc",      # Sea level pressure
    :relative_humidity               => "RYF.rhuss.1990_1991.nc",    # Surface relative humidity
    :downwelling_longwave_radiation  => "RYF.rlds.1990_1991.nc",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "RYF.rsds.1990_1991.nc",     # Downwelling shortwave radiation
    :temperature                     => "RYF.tas.1990_1991.nc",      # Near-surface air temperature
    :eastward_velocity               => "RYF.uas.1990_1991.nc",      # Eastward near-surface wind
    :northward_velocity              => "RYF.vas.1990_1991.nc",      # Northward near-surface wind
)

jra55_short_names = Dict(
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

field_time_series_short_names = Dict(
    :river_freshwater_flux           => "Fri", # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "Fra", # Freshwater flux from rainfall
    :snow_freshwater_flux            => "Fsn", # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "Fic", # Freshwater flux from calving icebergs
    :specific_humidity               => "qa",  # Surface specific humidity
    :sea_level_pressure              => "pa",  # Sea level pressure
    :relative_humidity               => "rh",  # Surface relative humidity
    :downwelling_longwave_radiation  => "Ql",  # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "Qs",  # Downwelling shortwave radiation
    :temperature                     => "Ta",  # Near-surface air temperature
    :eastward_velocity               => "ua",  # Eastward near-surface wind
    :northward_velocity              => "va",  # Northward near-surface wind
)

urls = Dict(
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

# URLs for the JRA55 datasets specific to each version
function urls(metadata::JRA55Metadata)
    return "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/monthly/"
end

function download_dataset!(metadata::JRA55Metadata;
                           url = urls(metadata))

    for data in metadata
        filename  = metadata_filename(data)
        shortname = short_name(data)

        if !isfile(filename)
            fileurl = joinpath(url, shortname, year, filename)
            download(url, filepath)
        end
    end

    return nothing
end