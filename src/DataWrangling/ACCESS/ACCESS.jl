module ACCESS

export ACCESS1deg, ACCESS025deg

using Downloads
using Oceananigans
using Scratch
using Oceananigans.DistributedComputations: @root

using ..DataWrangling: Metadatum, metadata_path

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


download_ACCESS_cache::String = ""
function __init__()
    global download_ACCESS_cache = @get_scratch!("ACCESS")
end
    

ACCESS_bathymetry_variable_names = Dict(
    :bottom_height => "depth",
)

struct ACCESS1deg end
struct ACCESS025deg end

const ACCESSBathymetry = Union{ACCESS1deg, ACCESS025deg}

default_download_directory(::ACCESSBathymetry) = download_ACCESS_cache

reversed_vertical_axis(::ACCESSBathymetry) = false

longitude_interfaces(::ACCESSBathymetry) = (-280, 80)
latitude_interfaces(::ACCESSBathymetry) = (-90, 90)

Base.size(::ACCESS1deg) = (300, 360, 1)
Base.size(::ACCESS025deg) = (1142, 1440, 1)
Base.size(dataset::ACCESSBathymetry, variable) = size(dataset)

all_dates(::ACCESSBathymetry, args...) = nothing
first_date(::ACCESSBathymetry, args...) = nothing
last_date(::ACCESSBathymetry, args...) = nothing

const ACCESSMetadatum = Metadatum{<:ACCESSBathymetry, <:Any, <:Any}

dataset_variable_name(data::ACCESSMetadatum) = ACCESS_bathymetry_variable_names[data.name]

z_interfaces(::ACCESSMetadatum) = (0, 1)
metadata_url(::ACCESSMetadatum) = nothing

metadata_filename(metadatum::Metadatum{ACCESS1deg, <:Any, <:Any}) = "access-om3-topo-1deg.nc"
metadata_filename(metadatum::Metadatum{ACCESS025deg, <:Any, <:Any}) = "access-om3-topo-0.25deg.nc"

function download_dataset(metadatum::ACCESSMetadatum)
    filepath = metadata_path(metadatum)

    @root if !isfile(filepath)
        throw(ArgumentError("$(metadata_filename(metadatum)) in $(filepath) not found!")) 
    end
    return filepath
end

end #module
