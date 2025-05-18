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
    # "single_column_os_papa_simulation.jl",
    # "one_degree_simulation.jl",
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
        # "Single-column ocean simulation" => "literated/single_column_os_papa_simulation.md",
        # "One-degree ocean simulation" => "literated/one_degree_simulation.md",
        # "Near-global ocean simulation" => "literated/near_global_ocean_simulation.md",
        ],

    "Interface fluxes" => "interface_fluxes.md",

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
