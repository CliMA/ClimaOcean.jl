module ECCO

export ECCOMetadatum, ECCO_immersed_grid, adjusted_ECCO_tracers, initialize!
export ECCO2Monthly, ECCO4Monthly, ECCO2Daily

using ClimaOcean
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: inpaint_mask!, NearestNeighborInpainting, download_progress,
                                compute_native_date_range

using Oceananigans
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

import ClimaOcean.DataWrangling: vertical_interfaces, variable_is_three_dimensional,
                                 shift_longitude_to_0_360, inpainted_metadata_path,
                                 longitude_shift

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

# ECCO4 data is on a -180, 180 longitude grid as opposed to ECCO2 data that
# is on a 0, 360 longitude grid. To make the data consistent, we shift ECCO4
# data by 180 degrees in longitude
longitude_shift(metadata::Metadata{<:ECCO4Monthly}) = 180

function inpainted_metadata_filename(metadata::ECCOMetadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::ECCOMetadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

end # Module
