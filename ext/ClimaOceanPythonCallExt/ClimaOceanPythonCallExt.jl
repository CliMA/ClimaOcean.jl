module ClimaOceanVerosExt

using ClimaOcean
using CondaPkg
using PythonCall
using Oceananigans
using Oceananigans.DistributedComputations: @root

using Dates: DateTime

include("VerosOceanSimulations/veros_ocean_simulation.jl")
include("veros_state_exchanger.jl")

end # module ClimaOceanVerosExt
