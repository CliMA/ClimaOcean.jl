module AtmosphereSimulations

using Oceananigans
using Oceananigans.Fields: Center
using Oceananigans.Grids: grid_name
using Oceananigans.OutputReaders: FieldTimeSeries, update_field_time_series!, extract_field_time_series
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Simulations: TimeStepWizard
using Oceananigans.Utils: prettysummary, Time

using Adapt
using Thermodynamics.Parameters: AbstractThermodynamicsParameters

import Oceananigans.TimeSteppers: time_step!, update_state!
import ClimaOcean: compute_net_atmosphere_fluxes!

include("atmospheric_thermodynamics_parameters.jl")
include("prescribed_atmosphere.jl")

end # module
