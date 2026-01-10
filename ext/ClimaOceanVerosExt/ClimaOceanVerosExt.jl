module ClimaOceanVerosExt

using ClimaOcean
using CondaPkg
using PythonCall
using Oceananigans

include("veros_ocean_simulation.jl")
include("veros_state_exchanger.jl")

end # module ClimaOceanVerosExt
