# ClimaOcean.jl

🌎 Realistic ocean-only and coupled ocean + sea-ice simulations driven by prescribed atmospheres and based on [Oceananigans](https://github.com/CliMA/Oceananigans.jl) and [ClimaSeaIce](https://github.com/CliMA/ClimaSeaIce.jl).

ClimaOcean implements a framework for driving coupled ocean and sea ice simulations with prescribed atmospheres, using bulk formula to compute atmosphere-ocean, atmosphere-ice, and ocean-ice fluxes of heat, momentum, and freshwater. ClimaOcean builds off the Oceananigans framework, which provides tools for gridded finite volume computations on CPUs and GPUs and building ocean-flavored fluid dynamics simulations. ClimaSeaIce, which provides software for both stand-alone and coupled sea ice simulations, is also built with Oceananigans.

ClimaOcean's core abstractions are `OceanSeaIceModel` which encapsulates the ocean simulation, sea ice simulation, prescribed atmospheric state, atmospheric thermodynamic parameters, and parameterizations that define how the three communicate. ClimaOcean also implements `ocean_simulation`, a utility for building realistic, hydrostatic ocean simulations with Oceananigans ensuring compatibility with `OceanSeaIceModel`.

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

!!! compat "Julia 1.9 is required"
    ClimaOcean requires Julia 1.9 or later.
    
## Quick start

The following script implements a near-global ocean simulation initialized from the [ECCO state estimate](https://gmd.copernicus.org/articles/8/3071/2015/) and coupled to a prescribed atmosphere derived from the [JRA55-do reanalysis](https://www.sciencedirect.com/science/article/pii/S146350031830235X):

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

# Build and run an OceanSeaIceModel (with no sea ice component) forced by JRA55 reanalysis
atmosphere = ClimaOcean.JRA55_prescribed_atmosphere(arch)
coupled_model = ClimaOcean.OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model, Δt=5minutes, stop_time=30days)
run!(simulation)
```

The simulation above achieves approximately 8 simulated years per day of wall time on an Nvidia H100 GPU.

We can leverage `Oceananigans` features to plot the surface speed at the end of the simulation:

```julia
u, v, w = ocean.model.velocities
speed = Field(sqrt(u^2 + v^2))
compute!(speed)

using GLMakie
heatmap(view(speed, :, :, ocean.model.grid.Nz), colorrange=(0, 0.5), colormap=:magma, nan_color=:lightgray)
```

![Uploading 379438791-4c484b93-38fe-4840-bf7d-63a3a59d29e1.png…]()
