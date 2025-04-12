using
  ClimaOcean,
  Documenter,
  Literate

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

to_be_literated = [
    "ecco_inspect_temperature_salinity.jl",
    # "generate_bathymetry.jl",
    # "generate_surface_fluxes.jl",
    # "single_column_os_papa_simulation.jl",
    # "one_degree_simulation.jl",
    # "mediterranean_simulation_with_ecco_restoring.jl",
    # "near_global_ocean_simulation.jl"
]

for file in to_be_literated
    filepath = joinpath(EXAMPLES_DIR, file)
    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor(), execute = true)
    end
end

#####
##### Build and deploy docs
#####

format = Documenter.HTML(collapselevel = 2,
                         size_threshold = nothing,
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         canonical = "https://clima.github.io/ClimaOceanDocumentation/dev/")

pages = [
    "Home" => "index.md",

    "Examples" => [
        # "Inspect ECCO2 data" => "literated/ecco_inspect_temperature_salinity.md",
        "Generate bathymetry" => "literated/generate_bathymetry.md",
        "Surface fluxes" => "literated/generate_surface_fluxes.md",
        "Single-column simulation" => "literated/single_column_os_papa_simulation.md",
        # "Mediterranean simulation with ECCO restoring" => "literated/mediterranean_simulation_with_ecco_restoring.md",
        "One-degree Ocean simulation" => "literated/one_degree_simulation.md",
        "Near-global Ocean simulation" => "literated/near_global_ocean_simulation.md",
        ],

    "Library" => [
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],
]

makedocs(sitename = "ClimaOcean.jl";
         format,
         pages,
         modules = [ClimaOcean],
         doctest = true,
         clean = true,
         warnonly = [:cross_references, :missing_docs],
         checkdocs = :exports)

@info "Clean up temporary .jld2 and .nc output created by doctests or literated examples..."

"""
    recursive_find(directory, pattern)

Return list of filepaths within `directory` that contains the `pattern::Regex`.
"""
recursive_find(directory, pattern) =
    mapreduce(vcat, walkdir(directory)) do (root, dirs, files)
        joinpath.(root, filter(contains(pattern), files))
    end

files = []
for pattern in [r"\.jld2", r"\.nc"]
    global files = vcat(files, recursive_find(@__DIR__, pattern))
end

for file in files
    rm(file)
end

ci_build = get(ENV, "CI", nothing) == "true"

if ci_build
    deploydocs(repo = "github.com/CliMA/ClimaOceanDocumentation.git",
               deploy_config = Documenter.Buildkite(),
               versions = ["stable" => "v^", "dev" => "dev", "v#.#.#"],
               forcepush = true,
               devbranch = "main",
               push_preview = true)
end
