using Test

@testset "ClimaOcean" begin
    include("test_configurations.jl")
    include("test_total_water_and_diagnostics.jl")
    include("test_one_degree_regression.jl")
end
