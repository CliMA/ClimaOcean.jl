module Atmospheres

export PrescribedAtmosphere, PrognosticAtmosphere

using KernelAbstractions: @kernel, @index

"""
    abstract type AbstractAtmosphere 

An abstract type representing an atmosphere configured to exchange fluxes with a sea ice model and 
and ocean model defined in ClimaOcean
"""
abstract type AbstractAtmosphere end

include("atmospheric_parameters.jl")
include("prescribed_atmospheres.jl")
include("prognostic_atmospheres.jl")


end