using ClimaOcean
using Oceananigans
using CUDA
using Test

using Oceananigans.Architectures: architecture
using Oceananigans.OutputReaders: interpolate!

using ClimaOcean

test_architectures = [CPU()]
