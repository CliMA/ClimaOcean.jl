pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add ClimaOcean to environment stack

using
  Documenter,
  Literate,
  ClimaOcean

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

to_be_literated = [
    "inspect_ecco_data.jl",
    # "generate_surface_fluxes.jl",
    # "single_column_simulation.jl",
    # "mediterranean_simulation_with_ecco_restoring.jl",
    # "near_global_omip_simulation.jl"
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

format = Documenter.HTML(
   collapselevel = 2,
  size_threshold = nothing,
      prettyurls = get(ENV, "CI", nothing) == "true",
       canonical = "https://clima.github.io/ClimaOceanDocumentation/dev/",
)

pages = [
    "Home" => "index.md",

    "Library" => [ 
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],

    "Examples" => [
        "Inspect ECCO2 data" => "literated/inspect_ecco_data.md",
        # "Surface fluxes" => "literated/generate_surface_fluxes.md",
        # "Single column simulation" => "literated/single_column_simulation.md",
        # "Mediterranean simulation with ECCO restoring" => "literated/mediterranean_simulation_with_ecco_restoring.md",
        # "Near global OMIP simulation" => "literated/near_global_omip_simulation.md",
        ]
]

makedocs(
     sitename = "ClimaOcean.jl",
      modules = [ClimaOcean],
       format = format,
        pages = pages,
      doctest = true,
        clean = true,
     warnonly = [:cross_references, :missing_docs],
    checkdocs = :exports
)

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

withenv("GITHUB_REPOSITORY" => "CliMA/ClimaOceanDocumentation") do
    documenter_home = ENV["HOME"]
    ENV["DOCUMENTER_KEY"] = readline("$HOME/.ssh/id_rsa")
    deploydocs(        repo = "github.com/CliMA/ClimaOceanDocumentation.git",
                   versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
                  forcepush = true,
                  devbranch = "main",
               push_preview = true)
end
