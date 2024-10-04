<!-- Title -->
<h1 align="center">
  ClimaOcean.jl
</h1>

<!-- description -->
<p align="center">
  <strong>🌎 Tools for building realistic ocean-only and coupled ocean + sea-ice simulations based on
          <a href="https://github.com/CliMA/Oceananigans.jl">Oceananigans</a>
          and <a href="https://github.com/CliMA/ClimaSeaIce.jl">ClimaSeaIce</a>.</strong>
</p>

<!-- Information badges -->
<p align="center">

   <a href="https://doi.org/10.5281/zenodo.7677442">
    <img alt="DOI" src="https://zenodo.org/badge/DOI/10.5281/zenodo.7677442.svg?style=flat-square">
  </a>

  <a href="https://clima.github.io/ClimaOceanDocumentation/dev">
    <img alt="Documentation" src="https://img.shields.io/badge/documentation-in%20development-orange?style=flat-square">
  </a>

  <a href="https://buildkite.com/clima/climaocean-ci">
    <img alt="Build status" src="https://badge.buildkite.com/3113cca353b83df3b5855d3f0d69827124614aef7017c835d2.svg?style=flat-square">
  </a>

</p>

## Installation

To install from a Julia REPL:

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/CliMA/ClimaOcean.jl.git")

julia> Pkg.instantiate()
```

Use `Pkg.add("url=https://github.com/CliMA/ClimaOcean.jl.git", rev="main")` to install the latest version of `ClimaOcean`.
For more information, see the [documentation for `Pkg.jl`](https://pkgdocs.julialang.org).

## Why? What's the difference between ClimaOcean and [Oceananigans](https://github.com/CliMA/Oceananigans.jl)?

[`Oceananigans`](https://github.com/CliMA/Oceananigans.jl) is a general-purpose package for ocean-flavored fluid dynamics. 
`ClimaOcean` _specializes_ `Oceananigans` for a specific application: realistic ocean simulations, and coupled ocean + sea ice simulations.

Realistic ocean simulations require

* Simulating the evolution of ocean temperature and salinity.
* Computing fluxes of heat, water vapor, momentum, and trace gases between the ocean and atmosphere, where the atmospheric state is either prescribed, or evolving in an atmosphere-ocean coupled configuration.
* Defining domain geometry using bathymetry observations.
* Initializing the ocean model with realistic temperature, salinity, and velocity fields.
    
`ClimaOcean` is built on top of `Oceananigans` and `ClimaSeaIce`, so it's important that `ClimaOcean` users become proficient with [`Oceananigans`](https://github.com/CliMA/Oceananigans.jl).
Note that `ClimaOcean` is currently focused on hydrostatic modeling with `Oceananigans`' `HydrostaticFreeSurfaceModel`.

    

