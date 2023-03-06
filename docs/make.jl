pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..")) # add ClimaOcean to environment stack

using
  Documenter,
  Literate,
  CairoMakie,  # so that Literate.jl does not capture precompilation output or warnings
  ClimaOcean

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

#####
##### Generate examples
#####

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/literated")

to_be_literated = []

for file in to_be_literated
    filepath = joinpath(EXAMPLES_DIR, file)
    Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor())
end

#####
##### Build and deploy docs
#####

format = Documenter.HTML(
  collapselevel = 2,
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
]

makedocs(
   sitename = "ClimaOcean.jl",
    modules = [ClimaOcean],
     format = format,
      pages = pages,
    doctest = true,
     strict = true,
      clean = true,
  checkdocs = :exports
)

@info "Cleaning up temporary .jld2 and .nc files created by doctests..."

for file in vcat(glob("docs/*.jld2"), glob("docs/*.nc"))
    rm(file)
end

withenv("GITHUB_REPOSITORY" => "CliMA/ClimaOceanDocumentation") do
    deploydocs(        repo = "github.com/CliMA/ClimaOceanDocumentation.git",
                   versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
                  forcepush = true,
                  devbranch = "main",
               push_preview = true)
end
