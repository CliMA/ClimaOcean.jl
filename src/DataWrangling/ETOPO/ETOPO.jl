module ETOPO

using Scratch
export ETOPOBathy

import ClimaOcean.DataWrangling:
    metadata_filename,
    default_download_directory,
    metadata_path,
    Metadata,
    Metadatum

import ClimaOcean.Bathymetry:
    download_bathymetry

download_ETOPO_cache::String = ""
function __init__()
    global download_ETOPO_cache = @get_scratch!("ETOPO")
end

ETOPO_bathymetry_variable_names = Dict(
    :temperature => "temperature",
    :salinity    => "salinity"
)

struct ETOPOBathy end

default_download_directory(::ETOPOBathy) = download_ETOPO_cache
Base.size(::ETOPOBathy, variable) = (360, 173, 42)
reversed_vertical_axis(::ETOPOBathy) = true
longitude_interfaces(::ETOPOBathy) = (0.5, 360.5)
latitude_interfaces(::ETOPOBathy) = (-83.5, 89.5)

const ETOPOMetadatum = Metadatum{<:ETOPOBathy}
const ETOPO_url  = "https://www.dropbox.com/scl/fi/6pwalcuuzgtpanysn4h6f/" *
"ETOPO_2022_v1_60s_N90W180_surface.nc?rlkey=2t7890ruyk4nd5t5eov5768lt&st=yfxsy1lu&dl=0"

function ETOPOMetadatum(name;
    dir = download_ETOPO_cache)
    return Metadatum(name; dir, dataset=ETOPOBathy())
end

function metadata_url(m::ETOPOMetadatum)
    return ETOPO_url
end

function metadata_filename(metadatum::ETOPOMetadatum)
    return "ETOPO_2022_v1_60s_N90W180_surface.nc"
end

function download_bathymetry(metadatum::ETOPOMetadatum)

    filepath = metadata_path(metadatum)
    fileurl  = metadata_url(metadatum)

    #TODO: embed this into a @root macro; see failed attempts in https://github.com/CliMA/ClimaOcean.jl/pull/391
    if !isfile(filepath)
        Downloads.download(fileurl, filepath; downloader, progress=download_progress)
    end

    return nothing
end

end #module