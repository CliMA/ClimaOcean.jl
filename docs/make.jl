using ClimaOcean
using Documenter
using DocumenterCitations
using Literate

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

bib_filepath = joinpath(dirname(@__FILE__), "climaocean.bib")
bib = CitationBibliography(bib_filepath, style=:authoryear)

#####
##### Example definition and filtering
#####

struct Example
    title::String
    basename::String
    full_year::Bool
end

const EXAMPLES_DIR   = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR     = joinpath(@__DIR__, "src/literated")

mkpath(OUTPUT_DIR)

# Examples from examples/ directory.
# Set `full_year = true` to run the full 2-year simulation.
# Set `full_year = false` to run only 100 time steps (for quick CI).
# Setting `CLIMAOCEAN_BUILD_ALL_EXAMPLES=true` overrides all examples to full year.
coupled_examples = [
    Example("Latitude-longitude",   "latitude_longitude_ocean_sea_ice",    false),
    Example("Half-degree tripolar", "half_degree_tripolar_ocean_sea_ice",  false),
    Example("One-degree tripolar",  "one_degree_tripolar_ocean_sea_ice",   false),
    Example("ORCA",                 "orca_ocean_sea_ice",                  false),
]

# # The 1/6° distributed simulation is run via MPI in a separate CI step.
# # This visualization example loads the saved output and is executed by Literate.
# distributed_examples = [
#     Example("Sixth-degree distributed", "visualize_sixth_degree_simulation", false),
# ]

# When CLIMAOCEAN_BUILD_ALL_EXAMPLES is set, override all examples to full year
build_all = get(ENV, "CLIMAOCEAN_BUILD_ALL_EXAMPLES", "false") == "true"

#####
##### Generate examples using Literate (each in a subprocess for memory isolation)
#####

# for example in vcat(coupled_examples, distributed_examples)
#     script_path = joinpath(EXAMPLES_DIR, example.basename * ".jl")
#     full_simulation = example.full_year || build_all
#     withenv("CLIMAOCEAN_FULL_SIMULATION" => string(full_simulation)) do
#         run(`$(Base.julia_cmd()) --color=yes --project=$(dirname(Base.active_project())) $(joinpath(@__DIR__, "literate.jl")) $(script_path) $(OUTPUT_DIR)`)
#     end
# end

#####
##### Build and deploy docs
#####

coupled_pages       = [ex.title => joinpath("literated", ex.basename * ".md") for ex in coupled_examples]
# distributed_pages   = [ex.title => joinpath("literated", ex.basename * ".md") for ex in distributed_examples]

format = Documenter.HTML(collapselevel = 2,
                         size_threshold = nothing,
                         canonical = "https://clima.github.io/ClimaOceanDocumentation/stable/")

pages = [
    "Home" => "index.md",

    "Ocean--sea ice simulations"  => coupled_pages,
    # "Distributed simulations"     => distributed_pages,

    "Library" => [
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
    ],

    "References" => "references.md",
]

modules = Module[]

for m in [ClimaOcean]
    if !isnothing(m)
        push!(modules, m)
    end
end

makedocs(; sitename = "ClimaOcean.jl",
         format, pages, modules,
         plugins = [bib],
         doctest = true,
         doctestfilters = [
             r"┌ Warning:.*",  # remove standard warning lines
             r"│ Use at own risk",
             r"└ @ .*",        # remove the source location of warnings
         ],
         clean = true,
         warnonly = [:cross_references, :missing_docs],
         checkdocs = :exports)

@info "Clean up temporary .jld2, .nc, and .mp4 output created by doctests or literated examples..."

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
