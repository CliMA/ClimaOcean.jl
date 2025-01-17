module Atmospheres

export PrescribedAtmosphere

using Oceananigans
using Oceananigans.Utils: KernelParameters
using Oceananigans.Grids: grid_name, architecture
using Oceananigans.Utils: prettysummary, launch!
using Oceananigans.Fields: Center
using Oceananigans.OutputReaders: FieldTimeSeries, update_field_time_series!, extract_field_time_series

using KernelAbstractions: @kernel, @index

include("atmospheric_parameters.jl")
include("prescribed_atmospheres.jl")

# Need to be extended by atmospheric models
surface_layer_height(atmos::PrescribedAtmosphere)      = atmos.reference_height
boundary_layer_height(atmos::PrescribedAtmosphere)     = atmos.boundary_layer_height
thermodynamics_parameters(atmos::PrescribedAtmosphere) = atmos.thermodynamics_parameters

end