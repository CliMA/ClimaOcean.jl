using ClimaOcean
using Oceananigans
using CUDA
using Test

using ClimaOcean.Bathymetry: download_bathymetry_cache
using ClimaOcean.DataWrangling
using ClimaOcean.EN4
using ClimaOcean.ECCO
using ClimaOcean.JRA55

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.OutputReaders: interpolate!

using CFTime
using Dates

using CUDA: @allowscalar

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]

start_date = DateTimeProlepticGregorian(1993, 1, 1)

# Fictitious grid that triggers bathymetry download
function download_bathymetry(; dir = download_bathymetry_cache,
                             filename = "ETOPO_2022_v1_60s_N90W180_surface.nc")

    grid = LatitudeLongitudeGrid(size = (10, 10, 1),
                                 longitude = (0, 100),
                                 latitude = (0, 50),
                                 z = (-6000, 0))

    bottom = regrid_bathymetry(grid; dir, filename)

    return nothing
end

test_datasets = (ECCO2Monthly(), ECCO2Daily(), ECCO4Monthly(), EN4Monthly())

test_ecco2_datasets = tuple((ds for ds in test_datasets if startswith(string(typeof(ds)), "ECCO2"))...)
test_ecco4_en4_datasets = tuple((ds for ds in test_datasets if !startswith(string(typeof(ds)), "ECCO2"))...)

test_ecco_datasets = tuple((ds for ds in test_datasets if startswith(string(typeof(ds)), "ECCO"))...)
