module ClimaOceanPythonCallExt

using ClimaOcean
using CondaPkg
using PythonCall
using Oceananigans
using Oceananigans.DistributedComputations: @root

using Dates: DateTime

include("clima_ocean_copernicus.jl")
include("clima_ocean_veros.jl")

end # module ClimaOceanPythonCallExt
