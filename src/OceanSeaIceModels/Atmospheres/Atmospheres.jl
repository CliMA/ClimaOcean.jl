module Atmospheres

export PrescribedAtmosphere

using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: KernelParameters
using Oceananigans.Grids: grid_name
using Oceananigans.Utils: prettysummary
using Oceananigans.Fields: Center
using Oceananigans.OutputReaders: FieldTimeSeries, update_field_time_series!, extract_field_time_series

include("atmospheric_parameters.jl")
include("prescribed_atmospheres.jl")

# Need to be extended by atmospheric models
surface_layer_height(atmos::PrescribedAtmosphere)      = atmos.reference_height
boundary_layer_height(atmos::PrescribedAtmosphere)     = atmos.boundary_layer_height
thermodynamics_parameters(atmos::PrescribedAtmosphere) = atmos.thermodynamics_parameters

end