include("runtests_setup.jl")

using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.EN4
using ClimaOcean.DataWrangling: NearestNeighborInpainting, metadata_path, native_times, download_dataset

using Dates
using Oceananigans.Grids: topology
using Oceananigans.OutputReaders: time_indices
using Oceananigans.TimeSteppers: update_state!
using Oceananigans.Units

using CUDA: @allowscalar

# Inpaint only the first two cells inside the missing mask
inpainting = NearestNeighborInpainting(2)

test_datasets = (ECCO2Daily(),)
test_ecco_datasets = tuple((ds for ds in test_datasets if startswith(string(typeof(ds)), "ECCO"))...)

include("test_ecco_en4.jl")
