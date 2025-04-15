module ECCO

export ECCOMetadatum, ECCO_field, ECCO_mask, ECCO_immersed_grid, adjusted_ECCO_tracers, initialize!
export ECCO2Monthly, ECCO4Monthly, ECCO2Daily
export ECCOFieldTimeSeries, ECCORestoring

using ClimaOcean
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: inpaint_mask!, NearestNeighborInpainting, download_progress,
                                compute_native_date_range, dataset_field
using ClimaOcean.InitialConditions: three_dimensional_regrid!

using Oceananigans
using Oceananigans: location
using Oceananigans.Architectures: architecture, child_architecture
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: DistributedField, all_reduce, barrier!
using Oceananigans.Utils

using KernelAbstractions: @kernel, @index
using NCDatasets
using JLD2
using Downloads: download
using Dates
using Adapt
using Scratch

import ClimaOcean.DataWrangling: vertical_interfaces, empty_field, variable_is_three_dimensional,
                                 shift_longitude_to_0_360, inpainted_metadata_path, default_inpainting

download_ECCO_cache::String = ""
function __init__()
    global download_ECCO_cache = @get_scratch!("ECCO")
end

include("ECCO_metadata.jl")
include("ECCO_mask.jl")

const SomeECCODataset = Union{ECCO2Monthly, ECCO4Monthly, ECCO2Daily}

vertical_interfaces(metadata::Metadata{<:SomeECCODataset}) =
    [
    -6128.75,
    -5683.75,
    -5250.25,
    -4839.75,
    -4452.25,
    -4087.75,
    -3746.25,
    -3427.75,
    -3132.25,
    -2859.75,
    -2610.25,
    -2383.74,
    -2180.13,
    -1999.09,
    -1839.64,
    -1699.66,
    -1575.64,
    -1463.12,
    -1357.68,
    -1255.87,
    -1155.72,
    -1056.53,
    -958.45,
    -862.10,
    -768.43,
    -678.57,
    -593.72,
    -515.09,
    -443.70,
    -380.30,
    -325.30,
    -278.70,
    -240.09,
    -208.72,
    -183.57,
    -163.43,
    -147.11,
    -133.45,
    -121.51,
    -110.59,
    -100.20,
    -90.06,
    -80.01,
    -70.0,
    -60.0,
    -50.0,
    -40.0,
    -30.0,
    -20.0,
    -10.0,
        0.0,
    ]

empty_ECCO_field(variable_name::Symbol; kw...) = empty_field(Metadatum(variable_name, dataset=ECCO4Monthly()); kw...)

# Only temperature and salinity need a thorough inpainting because of stability,
# other variables can do with only a couple of passes. Sea ice variables
# cannot be inpainted because zeros in the data are physical, not missing values.
function default_inpainting(metadata::ECCOMetadata)
    if metadata.name in [:temperature, :salinity]
        return NearestNeighborInpainting(Inf)
    elseif metadata.name in [:sea_ice_fraction, :sea_ice_thickness]
        return nothing
    else
        return NearestNeighborInpainting(5)
    end
end

# ECCO4 data is on a -180, 180 longitude grid as opposed to ECCO2 data that
# is on a 0, 360 longitude grid. To make the data consistent, we shift ECCO4
# data by 180 degrees in longitude
function shift_longitude_to_0_360(data, metadata::Metadata{<:ECCO4Monthly})
    Nx = size(data, 1)
    if variable_is_three_dimensional(metadata)
        shift = (Nx ÷ 2, 0, 0)
    else
        shift = (Nx ÷ 2, 0)
    end
    data = circshift(data, shift)

    return data
end

# Fallback
ECCO_field(var_name::Symbol; kw...) = ECCO_field(ECCOMetadata(var_name); kw...)

function inpainted_metadata_filename(metadata::ECCOMetadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::ECCOMetadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

include("ECCO_restoring.jl")

end # Module
