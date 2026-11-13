using ClimaOcean
using Documenter
using DocumenterCitations

bibliography = CitationBibliography(joinpath(@__DIR__, "src", "ClimaOcean.bib"))

makedocs(sitename = "ClimaOcean.jl",
         modules = [ClimaOcean],
         plugins = [bibliography],
         format = Documenter.HTML(collapselevel = 1, size_threshold = nothing),
         pages = ["Home" => "index.md",
                  "Library" => "library.md",
                  "References" => "references.md"])

deploydocs(repo = "github.com/CliMA/ClimaOcean.jl", push_preview = true)
