# ClimaOcean.jl

🌎 Ready-to-use ocean and ocean–sea ice configurations of the CliMA ocean model, for coupling with [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl) and for OMIP simulations that verify the model's biases and prediction skill. Built on [Oceananigans](https://github.com/CliMA/Oceananigans.jl), [ClimaSeaIce](https://github.com/CliMA/ClimaSeaIce.jl), and [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl).

ClimaOcean packages the CliMA ocean and sea-ice setups behind a small set of configuration constructors. The generic coupling framework, the air–sea/air–ice flux computations, and the dataset-wrangling utilities now live in [NumericalEarth](https://github.com/NumericalEarth/NumericalEarth.jl); ClimaOcean depends on it and re-exports its functionality (`OceanSeaIceModel`, `ocean_simulation`, the `Metadata`/`ECCO`/`JRA55` tooling, …), so scripts that use those names keep working.

The package serves three purposes:

1. **Ocean and sea-ice configurations** — turnkey ocean and sea-ice simulations at standard resolutions and grids.
2. **Coupling with ClimaCoupler** — those configurations are the ocean and sea-ice components that ClimaCoupler couples to a prognostic atmosphere and land model.
3. **OMIP simulations** — global ocean–sea ice runs following the OMIP protocol, with built-in diagnostics for verifying the mean state, biases, and prediction skill of the CliMA ocean model.

ClimaOcean is written in Julia by the [Climate Modeling Alliance](https://clima.caltech.edu)
and heroic external collaborators.

## Installation

ClimaOcean is a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/).

2. Launch Julia and type

```julia
julia> using Pkg

julia> Pkg.add("ClimaOcean")
```

!!! compat "Julia 1.10 is required"
    ClimaOcean requires Julia 1.10 or later.

## Ocean and sea-ice configurations

Each configuration returns an Oceananigans `Simulation` with realistic bathymetry, advection, closures, and air–sea flux boundary conditions already assembled. Building a coupled ocean + sea-ice setup is two lines:

```julia
using ClimaOcean

arch    = GPU()
ocean   = one_degree_tripolar_ocean(arch)
sea_ice = one_degree_tripolar_sea_ice(ocean)
```

The ocean configurations are [`latitude_longitude_ocean`](@ref), [`one_degree_tripolar_ocean`](@ref), [`half_degree_tripolar_ocean`](@ref), [`sixth_degree_tripolar_ocean`](@ref), [`tenth_degree_tripolar_ocean`](@ref), and [`orca_ocean`](@ref) (the NEMO eORCA mesh). The matching sea-ice configurations ([`one_degree_tripolar_sea_ice`](@ref), …, [`orca_sea_ice`](@ref)) build a prognostic `ClimaSeaIce` simulation on the ocean's grid. For memory-limited testing, [`simplified_ocean_closure`](@ref) replaces the full CATKE + Gent-McWilliams + biharmonic closure with a lightweight one.

Since `ocean.model` is an `Oceananigans.HydrostaticFreeSurfaceModel` and `sea_ice.model` is a `ClimaSeaIce.SeaIceModel`, the full Oceananigans/ClimaSeaIce toolset (initial conditions from `ECCO`/`EN4`, output writers, diagnostics) is available on the returned objects.

## Coupling with ClimaCoupler

The ocean and sea-ice configurations are the components that [ClimaCoupler](https://github.com/CliMA/ClimaCoupler.jl) couples to a prognostic atmosphere (e.g. ClimaAtmos) and land model. ClimaCoupler supplies the atmospheric state and drives the coupling clock, while ClimaOcean provides the ocean and sea ice together with the ocean ↔ sea-ice flux exchange (NumericalEarth's `compute_sea_ice_ocean_fluxes!`).

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

See the *Ocean–sea ice simulations* examples for complete, runnable scripts.

## OMIP simulations

[`omip_simulation`](@ref) builds a turnkey global ocean–sea ice run following the OMIP protocol (Griffies et al. 2016): a tripolar (or ORCA) grid forced by JRA55 reanalysis with salinity restoring, sea-ice initial conditions, and OMIP-protocol diagnostics already attached. These simulations are how the CliMA ocean model's mean state, biases, and prediction skill are verified.

```julia
using ClimaOcean
using Oceananigans.Units

simulation = omip_simulation(:halfdegree; arch = GPU(), Δt = 30minutes, stop_time = 2 * 365days)
run!(simulation)
```

The available configurations are `:halfdegree`, `:tenthdegree`, and `:orca`. ClimaOcean ships the diagnostics used to quantify drift and biases — [`add_omip_diagnostics!`](@ref) (SST/SSS/SSH/MLD, surface fluxes, sea-ice extent, zonal means, checkpoints) and [`strait_transports`](@ref) — together with a suite of interchangeable vertical-mixing closures for diagnostic comparisons (CATKE, KPP, NEMO-TKE, NORi, and Oceananigans' `RiBasedVerticalDiffusivity`).
