module Copernicus

using ClimaOcean.DataWrangling: Metadata, Metadatum
using Dates: DateTime, Day, Month

import Oceananigans.Fields:
    location

import ClimaOcean.DataWrangling:
    all_dates,
    default_download_directory,
    latitude_bounds,
    metadata_filename

using Scratch

download_Copernicus_cache::String = ""

function __init__()
    global download_Copernicus_cache = @get_scratch!("Copernicus")
end

# Datasets
abstract type CopernicusDataset end

default_download_directory(::CopernicusDataset) = download_Copernicus_cache

struct GLORYSDaily <: CopernicusDataset end
struct GLORYSMonthly <: CopernicusDataset end

all_dates(::GLORYSDaily, var) = range(DateTime("1993-01-01"), stop=DateTime("2025-03-24"), step=Day(1))
all_dates(::GLORYSMonthly, var) = range(DateTime("1993-01-01"), stop=DateTime("2024-12-01"), step=Month(1))

copernicusmarine_dataset_id(::GLORYSDaily) = "cmems_mod_glo_phy_my_0.083deg_P1D-m"
copernicusmarine_dataset_id(::GLORYSMonthly) = "cmems_mod_glo_phy_my_0.083deg_P1M-m"
# :static  => "cmems_mod_glo_phy_my_0.083deg_static",

struct CMEMSHourlyAnalysis <: CopernicusDataset end
copernicusmarine_dataset_id(::CMEMSHourlyAnalysis) = "cmems_mod_glo_phy_anfc_0.083deg_PT1H-m"

CopernicusMetadata{D} = Metadata{<:CopernicusDataset, D}
CopernicusMetadatum = Metadatum{<:CopernicusDataset}

copernicus_short_names = Dict(
    :temperature => "thetao",
    :salinity => "so",
    :sea_ice_concentration => "siconc",
    :sea_ice_thickness => "sithick",
    :u_velocity=> "uo",
    :v_velocity=> "vo",
    :sea_ice_u_velocity => "usi",
    :sea_ice_v_velocity => "vsi",
    :free_surface => "zos",
)       

start_date_str(date::DateTime) = string(date)
end_date_str(date::DateTime) = string(date)
start_date_str(dates) = first(dates) |> string
end_date_str(dates) = last(dates) |> string

function metadata_filename(metadata::CopernicusMetadata)
    shortname = copernicus_short_names[metadata.name]
    dataname = string(metadata.dataset)
    start_date = start_date_str(metadata.dates)
    end_date = end_date_str(metadata.dates)
    return string(shortname, "_",
                  dataname, "_",
                  start_date, "_",
                  end_date, ".nc")
end

latitude_bounds(::Metadata{<:CopernicusDataset}) = (-80, 90)
location(::Metadata{<:CopernicusDataset}) = (Center, Center, Center)

end # module Copernicus 
