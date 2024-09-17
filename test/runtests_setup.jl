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

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]
JRA55_data_directory = gpu_test ? "GPU_JRA55_data" : "CPU_JRA55_data"
bathymetry_data_directory = gpu_test ? "GPU_Bathymetry_data" : "CPU_Bathymetry_data"
ECCO_data_directory = gpu_test ? "GPU_ECCO_data" : "CPU_ECCO_data"  