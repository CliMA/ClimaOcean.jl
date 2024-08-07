# ClimaOcean.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7677442.svg)](https://doi.org/10.5281/zenodo.7677442)

[![Docs Dev](https://img.shields.io/badge/documentation-in%20development-orange)](https://clima.github.io/ClimaOceanDocumentation/dev)

[![Build status](https://badge.buildkite.com/3113cca353b83df3b5855d3f0d69827124614aef7017c835d2.svg)](https://buildkite.com/clima/climaocean-ci)

Tools for building and running realistic ocean-only and coupled ocean + sea-ice simulations based on [Oceananigans](https://github.com/CliMA/Oceananigans.jl) and [ClimaSeaIce](https://github.com/CliMA/ClimaSeaIce.jl).

## Installation

To install from a Julia REPL:

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/CliMA/ClimaOcean.jl.git")

julia> Pkg.instantiate()
```

Use `Pkg.add("url=https://github.com/CliMA/ClimaOcean.jl.git", rev="main") to install the latest version of `ClimaOcean`.
For more information, see the [documentation for `Pkg.jl`](https://pkgdocs.julialang.org).
