module Atmospheres

export PrescribedAtmosphere

using Oceananigans
using Oceananigans.Utils: KernelParameters
using Oceananigans.Grids: grid_name, architecture
using Oceananigans.Utils: prettysummary, launch!
using Oceananigans.Fields: Center, interpolate
using Oceananigans.OutputReaders: FieldTimeSeries, update_field_time_series!, extract_field_time_series

using KernelAbstractions: @kernel, @index

#####
##### Utility for interpolating tuples of fields
#####

# Note: assumes loc = (c, c, nothing) (and the third location should
# not matter.)
@inline interp_atmos_time_series(J, X, time, grid, args...) =
    interpolate(X, time, J, (c, c, nothing), grid, args...)

@inline interp_atmos_time_series(ΣJ::NamedTuple, args...) =
    interp_atmos_time_series(values(ΣJ), args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...) +
    interp_atmos_time_series(ΣJ[3], args...)

@inline interp_atmos_time_series(ΣJ::Tuple{<:Any, <:Any, <:Any, <:Any}, args...) =
    interp_atmos_time_series(ΣJ[1], args...) +
    interp_atmos_time_series(ΣJ[2], args...) +
    interp_atmos_time_series(ΣJ[3], args...) +
    interp_atmos_time_series(ΣJ[4], args...)


include("atmospheric_parameters.jl")
include("prescribed_atmospheres.jl")

# Need to be extended by atmospheric models
surface_layer_height(atmos::PrescribedAtmosphere)      = atmos.reference_height
boundary_layer_height(atmos::PrescribedAtmosphere)     = atmos.boundary_layer_height
thermodynamics_parameters(atmos::PrescribedAtmosphere) = atmos.thermodynamics_parameters

end