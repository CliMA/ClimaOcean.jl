# ClimaOcean.jl

ClimaOcean runs the global ocean--sea ice simulations of the Ocean Model Intercomparison Project
(OMIP-2; [Griffies et al. (2016)](@cite griffies2016omip), [Tsujino et al. (2020)](@cite tsujino2020evaluation))
with the ocean, sea ice and coupling of [NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl),
at three resolutions on the NEMO eORCA meshes:

| Configuration     | Mesh     | Eddies                          | Time step                     |
|:------------------|:---------|:--------------------------------|:------------------------------|
| `OneDegree()`     | eORCA1   | Gent-McWilliams and Redi        | 90 minutes                    |
| `QuarterDegree()` | eORCA025 | permitted                       | 4, then 30, then 45 minutes   |
| `TwelfthDegree()` | eORCA12  | resolved                        | 6, then 15 minutes            |

Each configuration carries the defaults of the corresponding production run, so that

```julia
using ClimaOcean
using Oceananigans.Units

simulation = omip_simulation(QuarterDegree())
simulation.stop_time = 300 * 365days
run!(simulation)
```

runs the quarter-degree OMIP simulation, forced by the multi-year JRA55-do reanalysis from 1958, and writes the OMIP
output of [`add_omip_diagnostics!`](@ref) together with checkpoints, from which `run!(simulation; pickup = true)`
restarts. A setting changes through the keyword of the same name, for example `OneDegree(biharmonic_timescale = 20days)`.

The scripts in `examples/` run the three configurations; the twelfth-degree one is distributed over four GPUs.
