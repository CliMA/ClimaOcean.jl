# Common test setup file to make stand-alone tests easy
using ClimaOcean
using Oceananigans
using CUDA
using Test

using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.OutputReaders: interpolate!
using Dates

using CUDA: @allowscalar

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]
start_date = DateTime(1993, 1, 1)

using CUDA

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

if test_group == :init || test_group == :all
end