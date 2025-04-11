<!-- Title -->
<h1 align="center">
  ClimaOcean.jl
</h1>

<!-- description -->
<p align="center">
  <strong>🌎 A framework for realistic ocean-only and coupled ocean + sea-ice simulations driven by prescribed atmospheres and based on <a href=https://github.com/CliMA/Oceananigans.jl>Oceananigans</a> and <a href=https://github.com/CliMA/ClimaSeaIce.jl>ClimaSeaIce</a></strong>.
</p>

###

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7677442.svg?style=flat-square)](https://doi.org/10.5281/zenodo.7677442)
[![Build status](https://badge.buildkite.com/3113cca353b83df3b5855d3f0d69827124614aef7017c835d2.svg?style=flat-square)](https://buildkite.com/clima/climaocean-ci)
[![Documentation](https://img.shields.io/badge/documentation-stable%20release-blue?style=flat-square)](https://clima.github.io/ClimaOceanDocumentation/stable)
[![Documentation](https://img.shields.io/badge/documentation-in%20development-orange?style=flat-square)](https://clima.github.io/ClimaOceanDocumentation/dev)

## Installation

ClimaOcean is a registered package. To install from a Julia REPL:

```julia
julia> using Pkg

julia> Pkg.add("ClimaOcean")

julia> Pkg.instantiate()
```

Use `Pkg.add(url="https://github.com/CliMA/ClimaOcean.jl.git", rev="main")` to install the latest version of `ClimaOcean`.
For more information, see the [documentation for `Pkg.jl`](https://pkgdocs.julialang.org).

## Why? What's the difference between ClimaOcean and [Oceananigans](https://github.com/CliMA/Oceananigans.jl)?

Oceananigans is a general-purpose library for ocean-flavored fluid dynamics.
ClimaOcean implements a framework for driving realistic Oceananigans simulations with prescribed atmospheres, and coupling them to prognostic sea ice simulations.

### A core abstraction: `ClimaOcean.OceanSeaIceModel`

Our system for realistic modeling is anchored by `ClimaOcean.OceanSeaIceModel`, which encapsulates the ocean simulation, sea ice simulation, prescribed atmospheric state, and specifies how the three communicate.
To illustrate how `OceanSeaIceModel` works we set up a simulation on a grid with 10 vertical levels and 1/4-degree horizontal resolution:

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
start_date = DateTime(1993, 1, 1)
set!(ocean.model,
     T=ClimaOcean.Metadata(:temperature; dates=start_date, dataset=ClimaOcean.ECCO4Monthly()),
     S=ClimaOcean.Metadata(:salinity;    dates=start_date, dataset=ClimaOcean.ECCO4Monthly()))

# Build and run an OceanSeaIceModel (with no sea ice component) forced by JRA55 reanalysis
atmosphere = ClimaOcean.JRA55PrescribedAtmosphere(arch)
coupled_model = ClimaOcean.OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model, Δt=5minutes, stop_time=30days)
run!(simulation)
```

The simulation above achieves approximately 8 simulated years per day of wall time on an Nvidia H100 GPU.

Since `ocean.model` is an `Oceananigans.HydrostaticFreeSurfaceModel`, we can leverage `Oceananigans` features in our scripts.
For example, to plot the surface speed at the end of the simulation we write

```julia
u, v, w = ocean.model.velocities
speed = Field(sqrt(u^2 + v^2))
compute!(speed)

using GLMakie
heatmap(view(speed, :, :, ocean.model.grid.Nz), colorrange=(0, 0.5), colormap=:magma, nan_color=:lightgray)
```

which produces

![image](https://github.com/user-attachments/assets/4c484b93-38fe-4840-bf7d-63a3a59d29e1)

### Additional features: a utility for `ocean_simulation`s and data wrangling

A second core abstraction in ClimaOcean is `ocean_simulation`. `ocean_simulation` configures an Oceananigans model for realistic simulations including temperature and salinity, the TEOS-10 equation of state, boundary conditions to store computed air-sea fluxes, the automatically-calibrated turbulence closure `CATKEVerticalDiffusivity`, and the [`WENOVectorInvariant` advection scheme](http://doi.org/10.1029/2023MS004130) for mesoscale-turbulence-resolving simulations.

ClimaOcean also provides convenience features for wrangling datasets of bathymetry, ocean temperature, salinity, ocean velocity fields, and prescribed atmospheric states.

ClimaOcean is built on top of Oceananigans and [ClimaSeaIce](https://github.com/CliMA/ClimaSeaIce.jl), so it's important that ClimaOcean users become proficient with [Oceananigans](https://github.com/CliMA/Oceananigans.jl).
Note that though ClimaOcean is currently focused on hydrostatic modeling with `Oceananigans.HydrostaticFreeSurfaceModel`, realistic nonhydrostatic modeling is also within the scope of this package.


### Citing

If you use ClimaOcean for your research, teaching, or fun 🤩, everyone in our community will be grateful
if you give credit by citing the corresponding Zenodo record, e.g.,

> Wagner, G. L. et al. (2025). CliMA/ClimaOcean.jl: v0.5.4 (v0.5.4). Zenodo. https://doi.org/10.5281/zenodo.15042648

and also the recent [preprint submitted to the Journal of Advances in Modeling Earth Systems](https://arxiv.org/abs/2502.14148) that presents an overview of all the things that make Oceananigans unique:

> "High-level, high-resolution ocean modeling at all scales with Oceananigans"
>
> by Gregory L. Wagner, Simone Silvestri, Navid C. Constantinou, Ali Ramadhan, Jean-Michel Campin,
> Chris Hill, Tomas Chor, Jago Strong-Wright, Xin Kai Lee, Francis Poulin, Andre Souza, Keaton J. Burns,
> John Marshall, Raffaele Ferrari
>
> submitted to the Journal of Advances in Modeling Earth Systems, arXiv 2502.14148

<details><summary>bibtex</summary>
  <pre><code>@article{Oceananigans-overview-paper-2025,
  title = {{High-level, high-resolution ocean modeling at all scales with Oceananigans}},
  author = {G. L. Wagner and S. Silvestri and N. C. Constantinou and A. Ramadhan and J.-M. Campin and C. Hill and T. Chor and J. Strong-Wright and X. K. Lee and F. Poulin and A. Souza and K. J. Burns and J. Marshall and R. Ferrari},
  journal = {arXiv preprint},
  year = {2025},
  archivePrefix = {arXiv},
  eprint = {2502.14148},
  doi = {10.48550/arXiv.2502.14148},
  notes = {submitted to the Journal of Advances in Modeling Earth Systems},
}</code></pre>
</details>
