module VerosOceanSimulations

using CondaPkg

using Oceananigans.Grids: topology
using ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity, SeaIceSimulation

import Oceananigans.Fields: set!
import Oceananigans.TimeSteppers: time_step!, initialize!

import ClimaOcean.OceanSeaIceModels: OceanSeaIceModel, default_nan_checker
import Oceananigans.Architectures: architecture

import Base: eltype

"""
    install_veros()

Install the Veros ocean model Marine CLI using CondaPkg.
Returns a NamedTuple containing package information if successful.
"""
function install_veros()
    CondaPkg.add_pip("veros")
    cli = CondaPkg.which("veros")
    @info "... the veros CLI has been installed at $(cli)."
    return cli
end

include("veros_ocean_simulation.jl")
include("veros_state_exchanger.jl")

end # module VerosOceanSimulations
