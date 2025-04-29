module Copernicus

using NCDatasets

using Oceananigans.Fields: Center
using ClimaOcean.DataWrangling: Metadata, Metadatum, metadata_path
using Dates: DateTime, Day, Month

import Oceananigans.Fields:
    location

import ClimaOcean.DataWrangling:
    all_dates,
    dataset_variable_name,
    default_download_directory,
    longitude_shift,
    absolute_latitude_bounds,
    metadata_filename,
    inpainted_metadata_path,
    reversed_vertical_axis,
    mangle_data,
    vertical_interfaces

using Scratch

download_Copernicus_cache::String = ""
function __init__()
    global download_Copernicus_cache = @get_scratch!("Copernicus")
end

# Datasets
abstract type CopernicusDataset end

default_download_directory(::CopernicusDataset) = download_Copernicus_cache

# This contains "static" variables -- eg the grid
struct GLORYSStatic <: CopernicusDataset end
struct GLORYSDaily <: CopernicusDataset end
struct GLORYSMonthly <: CopernicusDataset end

dataset_name(::GLORYSStatic) = "GLORYSStatic"
dataset_name(::GLORYSDaily) = "GLORYSDaily"
dataset_name(::GLORYSMonthly) = "GLORYSMonthly"

all_dates(::GLORYSStatic, var) = [nothing]
all_dates(::GLORYSDaily, var) = range(DateTime("1993-01-01"), stop=DateTime("2025-03-24"), step=Day(1))
all_dates(::GLORYSMonthly, var) = range(DateTime("1993-01-01"), stop=DateTime("2024-12-01"), step=Month(1))

copernicusmarine_dataset_id(::GLORYSStatic) = "cmems_mod_glo_phy_my_0.083deg_static"
copernicusmarine_dataset_id(::GLORYSDaily) = "cmems_mod_glo_phy_my_0.083deg_P1D-m"
copernicusmarine_dataset_id(::GLORYSMonthly) = "cmems_mod_glo_phy_my_0.083deg_P1M-m"
# :static  => "cmems_mod_glo_phy_my_0.083deg_static",

struct CMEMSHourlyAnalysis <: CopernicusDataset end
copernicusmarine_dataset_id(::CMEMSHourlyAnalysis) = "cmems_mod_glo_phy_anfc_0.083deg_PT1H-m"

CopernicusMetadata{D} = Metadata{<:CopernicusDataset, D}
CopernicusMetadatum = Metadatum{<:CopernicusDataset}

Base.size(::CopernicusMetadatum) = (4320, 2040, 50, 1)
reversed_vertical_axis(::CopernicusMetadata) = true

copernicus_dataset_variable_names = Dict(
    :temperature => "thetao",
    :depth => "deptho",
    :salinity => "so",
    :sea_ice_concentration => "siconc",
    :sea_ice_thickness => "sithick",
    :u_velocity=> "uo",
    :v_velocity=> "vo",
    :sea_ice_u_velocity => "usi",
    :sea_ice_v_velocity => "vsi",
    :free_surface => "zos",
)       

start_date_str(date) = string(date)
end_date_str(date) = string(date)
start_date_str(dates::Vector) = first(dates) |> string
end_date_str(dates::Vector) = last(dates) |> string

dataset_variable_name(metadata::CopernicusMetadata) = copernicus_dataset_variable_names[metadata.name]

function metadata_prefix(metadata::CopernicusMetadata)
    var = copernicus_dataset_variable_names[metadata.name]
    dataset = dataset_name(metadata.dataset)
    start_date = start_date_str(metadata.dates)
    end_date = end_date_str(metadata.dates)
    return string(var, "_",
                  dataset, "_",
                  start_date, "_",
                  end_date)
end

function metadata_filename(metadata::CopernicusMetadata)
    prefix = metadata_prefix(metadata)
    return string(prefix, ".nc")
end

function inpainted_metadata_path(metadata::CopernicusMetadata)
    prefix = metadata_prefix(metadata)
    return string(prefix, "_inpainted.jld2")
end

absolute_latitude_bounds(::Metadata{<:CopernicusDataset}) = (-80, 90)
location(::Metadata{<:CopernicusDataset}) = (Center, Center, Center)

function vertical_interfaces(metadata::Metadata{<:CopernicusDataset})
    path = metadata_path(metadata)
    ds = Dataset(path)
    zc = - reverse(ds["depth"][:])
    dz = zc[2] - zc[1]
    zf = zc[1:end-1] .+ zc[2:end]
    push!(zf, 0)
    pushfirst!(zf, zf[1] - dz)
    return zf
end

longitude_shift(::CopernicusMetadata) = 180

end # module Copernicus 
