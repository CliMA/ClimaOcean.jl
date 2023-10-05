# Experiments to ClimaOcean

This repository holds ClimaOcean "experiments".
An experiment has a few components:

1. Script that builds and runs an `Oceananigans.Simulation`. More specifically, a "script"
   is a file (such as `script.jl`) that can be run at the command line via

```bash
$ julia --project script.jl`,
```

or from the Julia REPL via

```julia-repl
julia> include("script.jl")
```

2. An "environment" - aa script that builds and runs a `Simulation`


