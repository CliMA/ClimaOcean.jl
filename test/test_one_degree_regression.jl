include("one_degree_regression_setup.jl")
using Test

@testset "One-degree OMIP simulation reproduces its reference after 100 steps" begin
    simulation = one_degree_regression_simulation(regression_architecture)
    run!(simulation)

    @test iteration(simulation) == 100

    architecture_name = string(nameof(typeof(regression_architecture)))
    summary = regression_summary(simulation)

    # Printed so that a reference can be taken from the log of a run
    @info "One-degree regression summary on $architecture_name"
    TOML.print(Dict(architecture_name => summary); sorted = true)

    reference = TOML.parsefile(regression_reference_path)[architecture_name]

    for (name, value) in reference
        @testset "$name" begin
            @test summary[name] ≈ value rtol = 1e-8
        end
    end
end
