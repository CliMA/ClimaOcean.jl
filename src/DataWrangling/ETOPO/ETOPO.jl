module ETOPO

using Scratch
export ETOPOBathy

import ClimaOcean.DataWrangling:
    metadata_filename,
    default_download_directory,
    metadata_path,
    Metadata,
    Metadatum,
    all_dates,
    download_dataset

using Oceananigans.DistributedComputations: @root

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

all_dates(dataset::ETOPOBathy) = nothing
all_dates(dataset::ETOPOBathy, variable) = nothing

const ETOPOMetadatum = Metadata{<:ETOPOBathy, <:Any, <:Any} #Metadatum{<:ETOPOBathy}
const ETOPO_url  = "https://www.dropbox.com/scl/fi/6pwalcuuzgtpanysn4h6f/" *
"ETOPO_2022_v1_60s_N90W180_surface.nc?rlkey=2t7890ruyk4nd5t5eov5768lt&st=yfxsy1lu&dl=0"

# function ETOPOMetadatum(name;
#     dir = download_ETOPO_cache)
#     return Metadatum(name; dir, dataset=ETOPOBathy())
# end

function metadata_url(m::ETOPOMetadatum)
    return ETOPO_url
end

function metadata_filename(metadatum::ETOPOMetadatum)
    return "ETOPO_2022_v1_60s_N90W180_surface.nc"
end

function download_dataset(metadata::ETOPOMetadatum)
    dir = metadata.dir

    # Create a temporary directory to store the .netrc file
    # The directory will be deleted after the download is complete
    @root mktempdir(dir) do tmp

        # Write down the username and password in a .netrc file
        ntasks = Threads.nthreads()

        asyncmap(metadata; ntasks) do metadatum # Distribute the download among tasks

            fileurl  = metadata_url(metadatum)
            filepath = metadata_path(metadatum)

            if !isfile(filepath)
                @info "Downloading ETOPO data: $(metadatum.name) in $(metadatum.dir)..."
                Downloads.download(fileurl, filepath; downloader, progress=download_progress)
            end
        end
    end

    return nothing
end


end #module