using ClimaOcean
using Oceananigans
using CUDA
using Test

using Oceananigans.Architectures: architecture

test_architectures = [CPU()]
