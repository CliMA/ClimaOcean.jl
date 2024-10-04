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

`Oceananigans` is a general-purpose package for ocean-flavored fluid dynamics. 
`ClimaOcean` _specializes_ `Oceananigans` for a specific application: realistic ocean simulations, and coupled ocean + sea ice simulations.

To do this, `ClimaOcean` implements two core abstractions:
* `ocean_simulation`, and
* `OceanSeaIceModel`.

To illustrate how `ClimaOcean` extends `Oceananigans`, consider this simple one-layer near-global model with 1/4 degree resolution,

```julia
using ClimaOcean
using Oceananigans
using Oceananigans.Units

grid = LatitudeLongitudeGrid(size=(1440, 560, 1), longitude=(0, 360), latitude=(-70, 70), z=(-3000, 0))
bathymetry = regrid_bathymetry(grid)
grid = ImmersedBoundaryGrid(grid, PartialCellBottom(bathymetry))

ocean = ocean_simulation(; grid) # build ocean component with default advection schemes and turbulence closures
date  = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T=ECCOMetadata(:temperature; date), S=ECCOMetadata(:salinity; date))

atmosphere = JRA55PrescribedAtmosphere() # a prescribed atmosphere based on JRA55 reanalysis
coupled_model = OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model, Δt=10minutes, stop_time=30days)
run!(simulation)
```

`ocean_simulation` configures an Oceananigans model for realistic simulations including temperature and salinity, the TEOS-10 equation of state, boundary conditions to store computed air-sea fluxes, the automatically-calibrated turbulence closure `CATKEVerticalDiffusivity`, and the [`WENOVectorInvariant` advection scheme](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023MS004130) for mesoscale-turbulence-resolving simulations.

`OceanSeaIceModel` provides a framework for coupled modeling that encapsulates the ocean simulation, a prescribed atmosphere, and optionally, a sea ice simulation.
`OceanSeaIceModel` is tasked with the computation of air-sea, air-ice, and ice-ocean fluxes.

In addition to these core abstractions `ClimaOcean` provides convenience features for wrangling datasets of bathymetry, ocean temperature, salinity, and velocity fields, and prescribed atmospheric states.
    
`ClimaOcean` is built on top of `Oceananigans` and `ClimaSeaIce`, so it's important that `ClimaOcean` users become proficient with [`Oceananigans`](https://github.com/CliMA/Oceananigans.jl).
Note that `ClimaOcean` is currently focused on hydrostatic modeling with `Oceananigans`' `HydrostaticFreeSurfaceModel`.

    

