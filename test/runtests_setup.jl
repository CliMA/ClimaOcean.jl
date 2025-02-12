using ClimaOcean
using Oceananigans
using CUDA
using Test

using ClimaOcean.DataWrangling
using ClimaOcean.ECCO
using ClimaOcean.JRA55

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.OutputReaders: interpolate!

using ClimaOcean
using ClimaOcean.Bathymetry: download_bathymetry_cache
using CFTime
using Dates 

using CUDA: @allowscalar

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]

# ECCO metadata for ECCO tests

start_date = DateTimeProlepticGregorian(1993, 1, 1)
end_date   = DateTimeProlepticGregorian(1993, 4, 1)
dates      = start_date : Month(1) : end_date

temperature_metadata = ECCOMetadata(:temperature, dates)
salinity_metadata    = ECCOMetadata(:salinity, dates)

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

