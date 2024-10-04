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

To illustrate how `ClimaOcean` extends `Oceananigans`, we set up a 10-layer near-global model at 1/4 degree resolution,

```julia
using Oceananigans
using Oceananigans.Units
using Dates, CFTime
import ClimaOcean

arch = GPU()
grid = LatitudeLongitudeGrid(arch,
                             size = (1440, 560, 10),
                             halo = (7, 7, 7),
                             longitude = (0, 360),
                             latitude = (-70, 70),
                             z = (-3000, 0))

bathymetry = ClimaOcean.regrid_bathymetry(grid) # builds gridded bathymetry based on ETOPO1
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))

# Build an ocean simulation initialized to the ECCO state estimate on Jan 1, 1993
ocean = ClimaOcean.ocean_simulation(grid)
date  = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T = ClimaOcean.ECCOMetadata(:temperature; date),
                  S = ClimaOcean.ECCOMetadata(:salinity; date))

# Build and run an OceanSeaIceModel (with no sea ice compoent this time) forced by JRA55 reanalysis
atmosphere = ClimaOcean.JRA55_prescribed_atmosphere(arch)
coupled_model = ClimaOcean.OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model, Δt=10minutes, stop_time=30days)
run!(simulation)
```

`ocean_simulation` configures an Oceananigans model for realistic simulations including temperature and salinity, the TEOS-10 equation of state, boundary conditions to store computed air-sea fluxes, the automatically-calibrated turbulence closure `CATKEVerticalDiffusivity`, and the [`WENOVectorInvariant` advection scheme](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023MS004130) for mesoscale-turbulence-resolving simulations.

`OceanSeaIceModel` provides a framework for coupled modeling that encapsulates the ocean simulation, a prescribed atmosphere, and optionally, a sea ice simulation.
`OceanSeaIceModel` is tasked with the computation of air-sea, air-ice, and ice-ocean fluxes.

In addition to these core abstractions `ClimaOcean` provides convenience features for wrangling datasets of bathymetry, ocean temperature, salinity, and velocity fields, and prescribed atmospheric states.
    
`ClimaOcean` is built on top of `Oceananigans` and `ClimaSeaIce`, so it's important that `ClimaOcean` users become proficient with [`Oceananigans`](https://github.com/CliMA/Oceananigans.jl).
Note that `ClimaOcean` is currently focused on hydrostatic modeling with `Oceananigans`' `HydrostaticFreeSurfaceModel`.

    

