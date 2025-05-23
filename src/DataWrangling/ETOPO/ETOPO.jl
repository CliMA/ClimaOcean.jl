module ETOPO

export ETOPO2022

using Downloads
using Oceananigans
using Oceananigans.DistributedComputations: @root
using Scratch

using ..DataWrangling: download_progress, Metadatum, metadata_path

import ClimaOcean.DataWrangling:
    metadata_filename,
    default_download_directory,
    all_dates,
    first_date,
    last_date,
    dataset_variable_name,
    download_dataset,
    longitude_interfaces,
    latitude_interfaces,
    z_interfaces,
    reversed_vertical_axis

download_ETOPO_cache::String = ""
function __init__()
    global download_ETOPO_cache = @get_scratch!("ETOPO")
end

ETOPO_bathymetry_variable_names = Dict(
    :bottom_height => "z",
)

struct ETOPO2022 end

default_download_directory(::ETOPO2022) = download_ETOPO_cache
reversed_vertical_axis(::ETOPO2022) = true
longitude_interfaces(::ETOPO2022) = (-180, 180)
latitude_interfaces(::ETOPO2022) = (-90, 90)
Base.size(::ETOPO2022) = (21600, 10800, 1)
Base.size(dataset::ETOPO2022, variable) = size(dataset)

all_dates(::ETOPO2022, args...) = nothing
first_date(::ETOPO2022, args...) = nothing
last_date(::ETOPO2022, args...) = nothing

const ETOPOMetadatum = Metadatum{<:ETOPO2022, <:Any, <:Any}

dataset_variable_name(data::ETOPOMetadatum) = ETOPO_bathymetry_variable_names[data.name]

const ETOPO_url = "https://www.dropbox.com/scl/fi/6pwalcuuzgtpanysn4h6f/" *
    "ETOPO_2022_v1_60s_N90W180_surface.nc?rlkey=2t7890ruyk4nd5t5eov5768lt&st=yfxsy1lu&dl=0"

z_interfaces(::ETOPOMetadatum) = (0, 1)
metadata_url(::ETOPOMetadatum) = ETOPO_url
metadata_filename(metadatum::ETOPOMetadatum) = "ETOPO_2022_v1_60s_N90W180_surface.nc"

function download_dataset(metadatum::ETOPOMetadatum)
    fileurl  = metadata_url(metadatum)
    filepath = metadata_path(metadatum)

    @root if !isfile(filepath)
        @info "Downloading ETOPO data: $(metadatum.name) in $(metadatum.dir)..."
        Downloads.download(fileurl, filepath; progress=download_progress)
    end
    return filepath
end

end #module
