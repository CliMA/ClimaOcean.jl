using ClimaOcean
using Oceananigans
using CUDA
using Test

using ClimaOcean.DataWrangling
using ClimaOcean.ECCO
using ClimaOcean.JRA55

using Oceananigans.Architectures: architecture
using Oceananigans.OutputReaders: interpolate!

test_architectures = [CPU()]
