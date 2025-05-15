module ETOPO

export ETOPOBathymetry

using Downloads
using Oceananigans
using Oceananigans.DistributedComputations: @root
using Scratch

import ClimaOcean.DataWrangling:
    metadata_filename,
    default_download_directory,
    metadata_path,
    Metadata,
    Metadatum,
    all_dates,
    first_date,
    last_date,
    download_dataset,
    download_progress,
    longitude_interfaces,
    latitude_interfaces,
    is_three_dimensional,
    z_interfaces,
    reversed_vertical_axis

download_ETOPO_cache::String = ""
function __init__()
    global download_ETOPO_cache = @get_scratch!("ETOPO")
end

ETOPO_bathymetry_variable_names = Dict(
    :temperature => "temperature",
    :salinity    => "salinity"
)

struct ETOPOBathymetry end

default_download_directory(::ETOPOBathymetry) = download_ETOPO_cache
reversed_vertical_axis(::ETOPOBathymetry) = true
longitude_interfaces(::ETOPOBathymetry) = (-180, 180)
latitude_interfaces(::ETOPOBathymetry) = (-90, 90)
Base.size(::ETOPOBathymetry, variable) = (21600, 10800, 1)
Base.size(::ETOPOBathymetry) = (21600, 10800, 1)

all_dates(::ETOPOBathymetry) = nothing
all_dates(::ETOPOBathymetry, variable) = nothing
first_date(::ETOPOBathymetry, variable) = nothing
last_date(::ETOPOBathymetry, variable) = nothing

const ETOPOMetadatum = Metadatum{<:ETOPOBathymetry, <:Any, <:Any}

const ETOPO_url = "https://www.dropbox.com/scl/fi/6pwalcuuzgtpanysn4h6f/" *
    "ETOPO_2022_v1_60s_N90W180_surface.nc?rlkey=2t7890ruyk4nd5t5eov5768lt&st=yfxsy1lu&dl=0"

z_interfaces(::ETOPOMetadatum)= (0, 1)
metadata_url(::ETOPOMetadatum) = ETOPO_url
metadata_filename(metadatum::ETOPOMetadatum) = "ETOPO_2022_v1_60s_N90W180_surface.nc"

function download_dataset(metadatum::ETOPOMetadatum)
    dir = metadatum.dir

    # Create a temporary directory to store the .netrc file
    # The directory will be deleted after the download is complete
    fileurl  = metadata_url(metadatum)
    filepath = metadata_path(metadatum)

    if !isfile(filepath)
        @info "Downloading ETOPO data: $(metadatum.name) in $(metadatum.dir)..."
        Downloads.download(fileurl, filepath; progress=download_progress)
    end
    return filepath
end

end #module
