pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add ClimaOcean to environment stack

using Documenter
using Literate
using ClimaOcean

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

#####
##### Build and deploy docs
#####

format = Documenter.HTML(collapselevel = 2,
                         prettyurls = get(ENV, "CI", nothing) == "true",
                         canonical = "https://clima.github.io/ClimaOceanDocumentation/dev/")

pages = [
    "Home" => "index.md",

    "Library" => [ 
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
        ],
]

makedocs(sitename = "ClimaOcean.jl"; format, pages,
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

deploydocs(repo = "github.com/CliMA/ClimaOceanDocumentation.git",
           versions = ["stable" => "v^", "dev" => "dev", "v#.#.#"],
           forcepush = true,
           devbranch = "main",
           push_preview = true)
