# Run the one-degree regression simulation and store its summary as the reference of this architecture:
#
#     julia --project=test test/generate_one_degree_regression_reference.jl            # CPU
#     GPU_TEST=true julia --project=test test/generate_one_degree_regression_reference.jl

include("one_degree_regression_setup.jl")

simulation = one_degree_regression_simulation(regression_architecture)
run!(simulation)

references = isfile(regression_reference_path) ? TOML.parsefile(regression_reference_path) : Dict{String, Any}()
references[string(nameof(typeof(regression_architecture)))] = regression_summary(simulation)

open(io -> TOML.print(io, references; sorted = true), regression_reference_path, "w")
