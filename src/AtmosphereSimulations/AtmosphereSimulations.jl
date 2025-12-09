module AtmosphereSimulations

using Oceananigans
using Oceananigans.Fields: Center
using Oceananigans.Grids: grid_name
using Oceananigans.OutputReaders: FieldTimeSeries, update_field_time_series!, extract_field_time_series
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Utils: prettysummary, Time

using Adapt
using Thermodynamics.Parameters: AbstractThermodynamicsParameters
using ClimaOcean.OceanSeaIceModels: OceanSeaIceModels

import Oceananigans.TimeSteppers: time_step!, update_state!


include("thermodynamics_parameters.jl")
include("prescribed_atmosphere.jl")
include("prescribed_atmosphere_exchanger.jl")
include("interpolate_atmospheric_state.jl")

end # module AtmosphereSimulations
