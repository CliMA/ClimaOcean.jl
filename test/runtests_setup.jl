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
end_date   = DateTimeProlepticGregorian(1993, 4, 1)
dates      = start_date : Month(1) : end_date

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

# Trigger downloading JRA55
arch = first(test_architectures)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))

test_ECCO_datasets = (ECCO4Monthly(), ECCO2Daily(), ECCO2Monthly())
test_datasets = (test_ECCO_datasets..., EN4Monthly())
