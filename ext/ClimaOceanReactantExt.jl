module ClimaOceanReactantExt

using Reactant
using Oceananigans.Architectures: ReactantState
using Oceananigans.DistributedComputations: Distributed

using ClimaOcean: OceanSeaIceModel

import Oceananigans

const OceananigansReactantExt = Base.get_extension(
     Oceananigans, :OceananigansReactantExt
)

const ReactantOSIM{I, A, O, F, C} = Union{
    OceanSeaIceModel{I, A, O, F, C, <:ReactantState},
    OceanSeaIceModel{I, A, O, F, C, <:Distributed{ReactantState}},
}

end # module ClimaOceanReactantExt
