# ClimaOcean.jl

🌎 Realistic ocean-only and coupled ocean + sea-ice simulations driven by prescribed atmospheres and based on [Oceananigans](https://github.com/CliMA/Oceananigans.jl) and [ClimaSeaIce](https://github.com/CliMA/ClimaSeaIce.jl).

ClimaOcean implements a framework for driving coupled ocean and sea ice simulations with prescribed atmospheres, using bulk formula to compute atmosphere-ocean, atmosphere-ice, and ocean-ice fluxes of heat, momentum, and freshwater. ClimaOcean builds off the Oceananigans framework, which provides tools for gridded finite volume computations on CPUs and GPUs and building ocean-flavored fluid dynamics simulations. ClimaSeaIce, which provides software for both stand-alone and coupled sea ice simulations, is also built with Oceananigans.

ClimaOcean's core abstractions are `OceanSeaIceModel` which encapsulates the ocean simulation, sea ice simulation, prescribed atmospheric state, atmospheric thermodynamic parameters, and parameterizations that define how the three communicate. ClimaOcean also implements `ocean_simulation`, a utility for building realistic, hydrostatic ocean simulations with Oceananigans ensuring compatibility with `OceanSeaIceModel`.

ClimaOcean is written in Julia by the [Climate Modeling Alliance](https://clima.caltech.edu)
and heroic external collaborators.

## Quick install

ClimaOcean is a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/).

2. Launch Julia and type

```julia
julia> using Pkg

julia> Pkg.add("ClimaOcean")
```

!!! compat "Julia 1.9 is required"
    ClimaOcean requires Julia 1.9 or later.
