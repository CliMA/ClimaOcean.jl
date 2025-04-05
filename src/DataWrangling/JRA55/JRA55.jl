module JRA55

export JRA55FieldTimeSeries, JRA55PrescribedAtmosphere, RepeatYearJRA55, MultiYearJRA55

using Oceananigans
using Oceananigans.Units
 
using Oceananigans.Architectures: arch_array
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: λnodes, φnodes, on_architecture
using Oceananigans.Fields: interpolate!
using Oceananigans.OutputReaders: Cyclical, TotallyInMemory, AbstractInMemoryBackend, FlavorOfFTS, time_indices

using ClimaOcean

using ClimaOcean.OceanSeaIceModels:
    PrescribedAtmosphere,
    TwoBandDownwellingRadiation

using CUDA: @allowscalar

using NCDatasets
using JLD2 
using Dates
using Scratch

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend, update_field_time_series!
using Downloads: download

download_JRA55_cache::String = ""

function __init__()
    global download_JRA55_cache = @get_scratch!("JRA55")
end

include("JRA55_metadata.jl")
include("JRA55_field_time_series.jl")
include("JRA55_prescribed_atmosphere.jl")

end # module
