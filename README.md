<!-- Title -->
<h1 align="center">
  ClimaOcean.jl
</h1>

<!-- description -->
<p align="center">
  <strong>🌎 Ready-to-use ocean and ocean + sea-ice configurations of the CliMA ocean model — for coupling with <a href=https://github.com/CliMA/ClimaCoupler.jl>ClimaCoupler</a> and for OMIP simulations that verify the model's biases and prediction skill. Built on <a href=https://github.com/CliMA/Oceananigans.jl>Oceananigans</a>, <a href=https://github.com/CliMA/ClimaSeaIce.jl>ClimaSeaIce</a>, and <a href=https://github.com/NumericalEarth/NumericalEarth.jl>NumericalEarth</a></strong>.
</p>

###

> [!IMPORTANT]
> The generic coupling framework and data-wrangling utilities originally developed in ClimaOcean now live in [**NumericalEarth.jl**](https://github.com/NumericalEarth/NumericalEarth.jl), a package for building coupled Earth system models with interchangeable components.
> ClimaOcean depends on NumericalEarth and re-exports its functionality (`OceanSeaIceModel`, `ocean_simulation`, the `Metadata`/`ECCO`/`JRA55` data tooling, etc.), so existing scripts that use these names keep working.
> ClimaOcean itself is now focused on **realistic ocean and ocean + sea-ice configurations** of the CliMA ocean model.
> See the [discussion](https://github.com/CliMA/ClimaOcean.jl/discussions/675) for more details.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7677442.svg?style=flat-square)](https://doi.org/10.5281/zenodo.7677442)
[![Build status](https://badge.buildkite.com/3113cca353b83df3b5855d3f0d69827124614aef7017c835d2.svg?style=flat-square)](https://buildkite.com/clima/climaocean-ci)
[![Documentation](https://img.shields.io/badge/documentation-stable%20release-blue?style=flat-square)](https://clima.github.io/ClimaOceanDocumentation/stable/)
[![Documentation](https://img.shields.io/badge/documentation-in%20development-orange?style=flat-square)](https://clima.github.io/ClimaOceanDocumentation/dev/)

## Installation

ClimaOcean is a registered package. To install from a Julia REPL:

```julia
julia> using Pkg

julia> Pkg.add("ClimaOcean")

julia> Pkg.instantiate()
```

Use `Pkg.add(url="https://github.com/CliMA/ClimaOcean.jl.git", rev="main")` to install the latest version of `ClimaOcean`.
For more information, see the [documentation for `Pkg.jl`](https://pkgdocs.julialang.org).

## What ClimaOcean provides

ClimaOcean packages the CliMA ocean and sea-ice setups behind three small interfaces.

### 1. Ocean and sea-ice configurations

Each configuration returns an `Oceananigans` `Simulation` with realistic bathymetry, advection, closures, and air–sea flux boundary conditions already assembled, so a coupled ocean + sea-ice setup is two lines:

```julia
using ClimaOcean

arch    = GPU()
ocean   = one_degree_tripolar_ocean(arch)
sea_ice = one_degree_tripolar_sea_ice(ocean)
```

Ocean configurations: `latitude_longitude_ocean`, `one_degree_tripolar_ocean`, `half_degree_tripolar_ocean`, `sixth_degree_tripolar_ocean`, `tenth_degree_tripolar_ocean`, and `orca_ocean` (NEMO eORCA mesh).
Matching sea-ice configurations (`one_degree_tripolar_sea_ice`, …, `orca_sea_ice`) build a prognostic `ClimaSeaIce` simulation on the ocean's grid.
For memory-limited testing, `simplified_ocean_closure()` swaps the full CATKE + Gent-McWilliams + biharmonic closure for a lightweight one.

Because `ocean.model` is an `Oceananigans.HydrostaticFreeSurfaceModel` and `sea_ice.model` is a `ClimaSeaIce.SeaIceModel`, the full Oceananigans/ClimaSeaIce toolset (initial conditions from `ECCO`/`EN4`, output writers, diagnostics) is available on the returned objects.

### 2. Coupling with ClimaCoupler

The ocean and sea-ice configurations are the components that [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl) couples to a prognostic atmosphere (e.g. ClimaAtmos) and land model.
ClimaCoupler supplies the atmospheric state and drives the coupling clock; ClimaOcean provides the ocean and sea ice and the ocean ↔ sea-ice flux exchange (via NumericalEarth's `compute_sea_ice_ocean_fluxes!`).
The same components can also be combined locally into a stand-alone coupled model driven by a prescribed atmosphere:

```julia
using ClimaOcean

arch       = GPU()
ocean      = orca_ocean(arch; closure = simplified_ocean_closure())
sea_ice    = orca_sea_ice(ocean)
atmosphere = JRA55PrescribedAtmosphere(arch)   # re-exported from NumericalEarth
radiation  = JRA55PrescribedRadiation(arch)

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
```

### 3. OMIP simulations

`omip_simulation` builds a turnkey global ocean–sea-ice run following the OMIP protocol (Griffies et al. 2016): a tripolar (or ORCA) grid forced by JRA55 reanalysis with salinity restoring, sea-ice initial conditions, and OMIP-protocol diagnostics already attached.
These simulations are how the CliMA ocean model's mean state, biases, and prediction skill are verified.

```julia
using ClimaOcean
using Oceananigans.Units

simulation = omip_simulation(:halfdegree; arch = GPU(), Δt = 30minutes, stop_time = 2 * 365days)
run!(simulation)
```

Configurations are `:halfdegree`, `:tenthdegree`, and `:orca`.
ClimaOcean ships the diagnostics used to quantify drift and biases — `add_omip_diagnostics!` (SST/SSS/SSH/MLD, surface fluxes, sea-ice extent, zonal means, checkpoints) and `strait_transports` — together with a suite of interchangeable vertical-mixing closures for diagnostic comparisons (CATKE, `KPP`, NEMO-TKE, NORi, and Oceananigans' `RiBasedVerticalDiffusivity`).

## Relationship to Oceananigans, ClimaSeaIce, and NumericalEarth

ClimaOcean is built on top of [Oceananigans](https://github.com/CliMA/Oceananigans.jl) (gridded finite-volume ocean dynamics on CPUs and GPUs), [ClimaSeaIce](https://github.com/CliMA/ClimaSeaIce.jl) (sea ice), and [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl) (the coupling framework, air–sea/air–ice flux computations, and dataset wrangling).
ClimaOcean users should become proficient with Oceananigans, since the returned models are Oceananigans models.
Though ClimaOcean currently focuses on hydrostatic modeling with `Oceananigans.HydrostaticFreeSurfaceModel`, realistic nonhydrostatic modeling is also within its scope.

## Citing

If you use ClimaOcean for your research, teaching, or fun 🤩, everyone in our community will be grateful
if you give credit by citing the corresponding Zenodo record, e.g.,

> Wagner, G. L. et al. (2025). CliMA/ClimaOcean.jl: v0.8.10 (v0.8.10). Zenodo. https://doi.org/10.5281/zenodo.7677442

and also the recent [preprint submitted to the Journal of Advances in Modeling Earth Systems](https://doi.org/10.48550/arXiv.2502.14148) that presents an overview of all the things that make Oceananigans unique:

> "High-level, high-resolution ocean modeling at all scales with Oceananigans"
>
> by Gregory L. Wagner, Simone Silvestri, Navid C. Constantinou, Ali Ramadhan, Jean-Michel Campin,
> Chris Hill, Tomas Chor, Jago Strong-Wright, Xin Kai Lee, Francis Poulin, Andre Souza, Keaton J. Burns,
> Siddhartha Bishnu, John Marshall, and Raffaele Ferrari
>
> submitted to the Journal of Advances in Modeling Earth Systems, arXiv:[2502.14148](https://doi.org/10.48550/arXiv.2502.14148)

<details><summary>bibtex</summary>
  <pre><code>@article{Oceananigans-overview-paper-2025,
  title = {{High-level, high-resolution ocean modeling at all scales with Oceananigans}},
  author = {G. L. Wagner and S. Silvestri and N. C. Constantinou and A. Ramadhan and J.-M. Campin and C. Hill and T. Chor and J. Strong-Wright and X. K. Lee and F. Poulin and A. Souza and K. J. Burns and S. Bishnu and J. Marshall and R. Ferrari},
  journal = {arXiv preprint},
  year = {2025},
  archivePrefix = {arXiv},
  eprint = {2502.14148},
  doi = {10.48550/arXiv.2502.14148},
  notes = {submitted to the Journal of Advances in Modeling Earth Systems},
}</code></pre>
</details>
