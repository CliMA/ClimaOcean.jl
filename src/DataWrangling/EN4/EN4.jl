module EN4

export EN4Metadatum, EN4_immersed_grid, adjusted_EN4_tracers, initialize!
export EN4Monthly

using ClimaOcean
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: inpaint_mask!, NearestNeighborInpainting, download_progress,
                                compute_native_date_range, Kelvin, Celsius

using Oceananigans
using Oceananigans.Architectures: architecture
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

import ClimaOcean.DataWrangling: vertical_interfaces, empty_field, is_three_dimensional,
                                 inpainted_metadata_path, longitude_shift

download_EN4_cache::String = ""
function __init__()
    global download_EN4_cache = @get_scratch!("EN4")
end

include("EN4_metadata.jl")

reversed_vertical_axis(::Metadata{<:EN4Monthly}) = true

vertical_interfaces(metadata::Metadata{<:EN4Monthly}) = [
    -5500.0,
    -5200.5986,
    -4901.459,
    -4602.695,
    -4304.469,
    -4007.0154,
    -3710.6636,
    -3415.8828,
    -3123.3313,
    -2833.926,
    -2548.9243,
    -2270.0168,
    -1999.4135,
    -1739.8875,
    -1494.7288,
    -1267.542,
    -1061.846,
    -880.5102,
    -725.1709,
    -595.8501,
    -490.9494,
    -407.6244,
    -342.3651,
    -291.5638,
    -251.9186,
    -220.6418,
    -195.5097,
    -174.8189,
    -157.3019,
    -142.0345,
    -128.3522,
    -115.7825,
    -103.9913,
    -92.7437,
    -81.8753,
    -71.2711,
    -60.8509,
    -50.5586,
    -40.3553,
    -30.214,
    -20.1158,
    -10.0475,
      0.0,
]

# EN4 data is shifted in longitude by 1 degree to the east.
# So to make the data consistent, we apply longitude shift.
longitude_shift(metadata::Metadata{<:EN4Monthly}) = 1

function inpainted_metadata_filename(metadata::EN4Metadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    var = string(metadata.name)
    return without_extension * "_" * var *"_inpainted.jld2"
end

inpainted_metadata_path(metadata::EN4Metadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

end # Module
