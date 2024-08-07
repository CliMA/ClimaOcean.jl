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

`ClimaOcean` is for _realistic_ ocean-only and ocean + sea-ice simulations, in a region of the ocean ("regional") or covering the whole Earth.
[Oceananigans](https://github.com/CliMA/Oceananigans.jl) is a lower-level package for simulating the dynamics of ocean-flavored fluids that can be used for _both_ idealized problems and, given enough effort, realistic problems as well.
While "idealized" problems come in multifarious shapes and sizes, "realistic" problems tend to be more narrowly defined, and require

* Simulating the evolution of specific tracers: ocean temperature (or heat), salinity, and sometimes ocean biogeochemistry.
* Computing fluxes of heat, water vapor, momentum, and trace gases between the ocean and atmosphere (where the atmospheric state is either prescribed or "coupled" and itself evolving) -- and also between sea ice and the atmosphere, when a sea ice component is included.
* Initializing the ocean model with realistic initial conditions derived from observations of the ocean, and realistic bathymetry.
    
`ClimaOcean` leverages `Oceananigans` and `ClimaSeaIce` to build `OceanSeaIceModel`s capable of meeting these requirements to simulate the dynamics of specific regions of the Earth's ocean.
So if you're using `ClimaOcean`, it's a very good idea to become proficient in [`Oceananigans`](https://github.com/CliMA/Oceananigans.jl) as well.
Note also that, at least at the moment, `ClimaOcean` is focused on hydrostatic modeling with `Oceananigans`' `HydrostaticFreeSurfaceModel`.

In summary, if you're interested in realistic, hydrostatic regional or global simulations you may find `ClimaOcean` useful.
Otherwise, you can stick with [Oceananigans](https://github.com/CliMA/Oceananigans.jl).

    

