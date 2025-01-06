module Atmospheres

export PrescribedAtmosphere, PrognosticAtmosphere

using KernelAbstractions: @kernel, @index

include("atmospheric_parameters.jl")
include("prescribed_atmospheres.jl")

# Need to be extended by atmospheric models
surface_layer_height(atmos::PrescribedAtmosphere)      = atmos.reference_height
boundary_layer_height(atmos::PrescribedAtmosphere)     = atmos.boundary_layer_height
thermodynamics_parameters(atmos::PrescribedAtmosphere) = atmos.thermodynamics_parameters

end