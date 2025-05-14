module ETOPO

using Scratch
export ETOPOBathymetry

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
    download_progress

using Downloads

using Oceananigans
using Oceananigans.DistributedComputations: @root

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
Base.size(::ETOPOBathymetry, variable) = (360, 173, 42)
reversed_vertical_axis(::ETOPOBathymetry) = true
longitude_interfaces(::ETOPOBathymetry) = (0.5, 360.5)
latitude_interfaces(::ETOPOBathymetry) = (-83.5, 89.5)

all_dates(dataset::ETOPOBathymetry) = nothing
all_dates(dataset::ETOPOBathymetry, variable) = nothing
function first_date(dataset::ETOPOBathymetry, variable)
    return nothing
end
last_date(dataset::ETOPOBathymetry, variable) = nothing

const ETOPOMetadatum = Metadatum{<:ETOPOBathymetry, <:Any, <:Any}
@info ETOPOMetadatum
const ETOPO_url  = "https://www.dropbox.com/scl/fi/6pwalcuuzgtpanysn4h6f/" *
"ETOPO_2022_v1_60s_N90W180_surface.nc?rlkey=2t7890ruyk4nd5t5eov5768lt&st=yfxsy1lu&dl=0"

function metadata_url(m::ETOPOMetadatum)
    return ETOPO_url
end

function metadata_filename(metadatum::ETOPOMetadatum)
    return "ETOPO_2022_v1_60s_N90W180_surface.nc"
end

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
    return nothing
end


end #module