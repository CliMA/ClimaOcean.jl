module ClimaOceanPythonCallExt

using ClimaOcean
using CondaPkg
using PythonCall
using Oceananigans
using Oceananigans.DistributedComputations: @root

using Dates: DateTime

include("copernicus.jl")
include("veros_ocean_simulation.jl")
include("veros_state_exchanger.jl")

end # module ClimaOceanPythonCallExt
