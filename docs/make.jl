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

# Set up a timer to print a space ' ' every 240 seconds. This is to avoid CI
# timing out when building demanding Literate.jl examples.
Timer(t -> println(" "), 0, interval=240)

format = Documenter.HTML(
  collapselevel = 2,
     prettyurls = get(ENV, "CI", nothing) == "true",
)

pages = [
    "Home" => "index.md",
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

withenv("GITHUB_REPOSITORY" => "CliMA/ClimaOcean.jl") do
    deploydocs(repo = "github.com/CliMA/ClimaOcean.jl.git",
               versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
               forcepush = true,
               devbranch = "main",
               push_preview = true)
end
