using CFTime
using Dates
import Dates: year, month, day
import Oceananigans.Fields: set!
import Base

# Ecco field used to set model's initial conditions
struct ECCOMetadata{D} 
    name  :: Symbol
    dates :: D
end

unprocessed_ecco4_file_prefix = Dict(
    :temperature  => ("OCEAN_TEMPERATURE_SALINITY_day_mean_", "_ECCO_V4r4_latlon_0p50deg.nc"),
    :salinity     => ("OCEAN_TEMPERATURE_SALINITY_day_mean_", "_ECCO_V4r4_latlon_0p50deg.nc"),
)

# We always start from 1992
ECCOMetadata(name::Symbol) = ECCOMetadata(name, DateTimeProlepticGregorian(1992, 1, 1))

# Treat Metadata as an array
Base.getindex(metadata::ECCOMetadata, i::Int) = @inbounds ECCOMetadata(metadata.name, metadata.dates[i])
Base.length(metadata::ECCOMetadata)       = length(metadata.dates)
Base.size(metadata::ECCOMetadata)         = size(metadata.dates)
Base.eltype(metadata::ECCOMetadata)       = Base.eltype(metadata.dates)
Base.axes(metadata::ECCOMetadata)         = Base.axes(metadata.dates)
Base.first(metadata::ECCOMetadata)        = @inbounds ECCOMetadata(metadata.name, metadata.dates[1])
Base.last(metadata::ECCOMetadata)         = @inbounds ECCOMetadata(metadata.name, metadata.dates[end])
Base.iterate(metadata::ECCOMetadata, i=1) = (@inline; (i % UInt) - 1 < length(metadata) ? (@inbounds ECCOMetadata(metadata.name, metadata.dates[i]), i + 1) : nothing)

Base.axes(metadata::ECCOMetadata{<:AbstractCFDateTime})  = 1
Base.first(metadata::ECCOMetadata{<:AbstractCFDateTime}) = metadata
Base.last(metadata::ECCOMetadata{<:AbstractCFDateTime})  = metadata

function date_string(metadata::ECCOMetadata{<:AbstractCFDateTime})
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    daystr   = string(Dates.day(metadata.dates), pad=2) 
    return "$(yearstr)-$(monthstr)-$(daystr)"
end

function file_name(metadata::ECCOMetadata{<:AbstractCFDateTime})
    variable_name = metadata.name
    prefix, postfix = unprocessed_ecco4_file_prefix[variable_name]
    datestr = date_string(metadata)
    return prefix * datestr * postfix
end

# Convenience functions
short_name(data::ECCOMetadata)     = ecco4_short_names[data.name]
field_location(data::ECCOMetadata) = ecco4_location[data.name]

variable_is_three_dimensional(data::ECCOMetadata) = 
    data.name == :temperature || 
    data.name == :salinity || 
    data.name == :u_velocity ||
    data.name == :v_velocity


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