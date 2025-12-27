using ClimaOcean
using Distributed

# Make sure the atmospheric data for the examples is downloaded
JRA55PrescribedAtmosphere() # This command downloads the atmospheric data used by the examples

# Distributed.addprocs(2)

# @everywhere begin
    using ClimaOcean
    using CUDA
    using Documenter
    using DocumenterCitations
    using Literate

    ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

    bib_filepath = joinpath(dirname(@__FILE__), "climaocean.bib")
    bib = CitationBibliography(bib_filepath, style=:authoryear)

    #####
    ##### Generate examples
    #####

    const EXAMPLES_DIR   = joinpath(@__DIR__, "..", "examples")
    const OUTPUT_DIR     = joinpath(@__DIR__, "src/literated")
    const DEVELOPERS_DIR = joinpath(@__DIR__, "src/developers")

    examples_pages = [
        "Python ocean forced simulation" => "literated/python_ocean_forced_simulation.md",
        "Single-column ocean simulation" => "literated/single_column_os_papa_simulation.md",
        "One-degree ocean--sea ice simulation" => "literated/one_degree_simulation.md",
        "Near-global ocean simulation" => "literated/near_global_ocean_simulation.md",
        "Global climate simulation" => "literated/global_climate_simulation.md",
    ]

    to_be_literated = map(examples_pages) do (_, mdpath)
        replace(basename(mdpath), ".md" => ".jl")
    end
# end

# Distributed.pmap(1:length(to_be_literated)) do n
for n in 1:length(to_be_literated)
    # device = Distributed.myid()
    # @info "switching to device $(device)"
    # CUDA.device!(device) # Set the correct GPU, the used GPUs will be number 2 and 3
    file = to_be_literated[n]
    filepath = joinpath(EXAMPLES_DIR, file)
    withenv("JULIA_DEBUG" => "Literate") do
        Literate.markdown(filepath, OUTPUT_DIR; flavor = Literate.DocumenterFlavor(), execute = true)
    end
    GC.gc(true)
    CUDA.reclaim()
end

# Distributed.rmprocs()

withenv("JULIA_DEBUG" => "Literate") do
    Literate.markdown(joinpath(DEVELOPERS_DIR, "slab_ocean.jl"), OUTPUT_DIR; flavor = Literate.DocumenterFlavor(), execute = true)
end

#####
##### Build and deploy docs
#####

format = Documenter.HTML(collapselevel = 2,
                         size_threshold = nothing,
                         canonical = "https://clima.github.io/ClimaOceanDocumentation/stable/")

pages = [
    "Home" => "index.md",

    "Examples" => examples_pages,

    "Developers" => [
        "OceanSeaIceModel interface" => "literated/slab_ocean.md",
        ],

    "Vertical grids" => "vertical_grids.md",

    "Metadata" => [
        "Overview" => "Metadata/metadata_overview.md",
        "Supported variables" => "Metadata/supported_variables.md",
    ],
    "Interface fluxes" => "interface_fluxes.md",

    "Library" => [
        "Contents"       => "library/outline.md",
        "Public"         => "library/public.md",
        "Private"        => "library/internals.md",
        "Function index" => "library/function_index.md",
    ],

    "References" => "references.md",
]

modules = Module[]
ClimaOceanSpeedyWeatherExt = isdefined(Base, :get_extension) ? Base.get_extension(ClimaOcean, :ClimaOceanSpeedyWeatherExt) : ClimaOcean.ClimaOceanSpeedyWeatherExt

for m in [ClimaOcean, ClimaOceanSpeedyWeatherExt]
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
