using ClimaOcean
using Oceananigans
using CUDA
using Test

using Oceananigans.Architectures: architecture, on_architecture
using Dates

using CUDA: @allowscalar

gpu_test = parse(Bool, get(ENV, "GPU_TEST", "false"))
test_architectures = gpu_test ? [GPU()] : [CPU()]

test_group = get(ENV, "TEST_GROUP", :all)
test_group = Symbol(test_group)

if test_group == :unit || test_group == :all
    include("test_module.jl")
    include("test_ocean_configurations.jl")
    include("test_sea_ice_configurations.jl")
    include("test_omip_configurations.jl")
end
