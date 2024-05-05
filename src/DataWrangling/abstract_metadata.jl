using CFTime
using Dates
import Dates: year, month, day
import Oceananigans.Fields: set!

abstract type AbstractMetadata end

# Ecco field used to set model's initial conditions
struct ECCOMetadata{D} <: AbstractMetadata
    name  :: Symbol
    dates :: D
end

unprocessed_ecco4_file_prefix = Dict(
    :temperature  => ("OCEAN_TEMPERATURE_SALINITY_day_mean_", "_ECCO_V4r4_latlon_0p50deg.nc"),
    :salinity     => ("OCEAN_TEMPERATURE_SALINITY_day_mean_", "_ECCO_V4r4_latlon_0p50deg.nc"),
)

# We always start from 1992
ECCOMetadata(name::Symbol) = ECCOMetadata(name, DateTimeAllLeap(1992, 1, 1))
date(data::ECCOMetadata)   = DateTimeAllLeap(data.year, data.month, data.day)

function extract_time_indices(metadata::ECCOMetadata) 
    dates = metadata.dates
    days  = dates .- dates[1] 
    return ceil.(Int, days ./ Day(1) .+ 1)
end

function date_string(metadata::ECCOMetadata{<:DateTimeAllLeap})
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    daystr   = string(Dates.day(metadata.dates), pad=2) 
    return "$(yearstr)-$(monthstr)-$(daystr)"
end

function file_name(metadata::ECCOMetadata{<:DateTimeAllLeap})
    variable_name = metadata.name
    prefix, postfix = unprocessed_ecco4_file_prefix[variable_name]
    datestr = date_string(metadata)
    return prefix * datestr * postfix
end

file_name(metadata::ECCOMetadata) = "ecco_$(metadata.name).nc"

variable_is_three_dimensional(data::ECCOMetadata) = 
    data.name == :temperature || 
    data.name == :salinity || 
    data.name == :u_velocity ||
    data.name == :v_velocity

short_name(data::ECCOMetadata) = data.name
field_location(data::ECCOMetadata)   = ecco4_location[data.name]

ecco4_short_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "UVEL",
    :v_velocity            => "VVEL",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_area_fraction => "SIarea"
)

ecco4_location = Dict(
    :temperature           => (Center, Center, Center),
    :salinity              => (Center, Center, Center),
    :sea_ice_thickness     => (Center, Center, Nothing),
    :sea_ice_area_fraction => (Center, Center, Nothing),
    :u_velocity            => (Face,   Center, Center),
    :v_velocity            => (Center, Face,   Center),
)

ecco4_remote_folder = Dict(
    :temperature           => "ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4",
    :salinity              => "ECCO_L4_TEMP_SALINITY_05DEG_DAILY_V4R4",
)

struct JRA55Metadata{I} <: AbstractMetadata
    name  :: Symbol
    time_indices :: I
end

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

extract_time_indices(metadata::JRA55Metadata) = metadata.time_indices

file_name(metadata::JRA55Metadata)  = filenames[metadata.name]
short_name(metadata::JRA55Metadata) = jra55_short_names[metadata.name]
variable_is_three_dimensional(metadata::JRA55Metadata) = false
field_location(::JRA55Metadata)     = (Center, Center, Nothing)

function download_dataset!(metadata::JRA55Metadata)
    filename = file_name(metadata)
    url = urls[metadata.name]
    isfile(filename) || download(url, filename)
    return nothing
end

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
