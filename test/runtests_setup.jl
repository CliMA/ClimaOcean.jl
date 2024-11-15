using ClimaOcean
using Oceananigans
using CUDA
using Test

using ClimaOcean.DataWrangling
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.JRA55: JRA55_field_time_series

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.OutputReaders: interpolate!

using ClimaOcean
using CFTime
using Dates 

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]

# ECCO metadata for ECCO tests

start_date = DateTimeProlepticGregorian(1993, 1, 1)
end_date   = DateTimeProlepticGregorian(1993, 4, 1)
dates      = start_date : Month(1) : end_date

temperature_metadata = ECCOMetadata(:temperature, dates)
salinity_metadata    = ECCOMetadata(:salinity, dates)

# Fictitious grid that triggers bathymetry download
function download_bathymetry()
    grid = LatitudeLongitudeGrid(size = (10, 10, 1), 
                                 longitude = (0, 100), 
                                 latitude = (0, 50),
                                 z = (-6000, 0))

    bottom = regrid_bathymetry(grid)

    return nothing
end
