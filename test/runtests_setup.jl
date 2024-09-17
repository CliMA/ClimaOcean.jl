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
data_directory = gpu_test ? "GPU_data" : "CPU_data"
