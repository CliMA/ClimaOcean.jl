module Atmospheres

export PrescribedAtmosphere, PrognosticAtmosphere

"""
    abstract type AbstractAtmosphere 

An abstract type representing an atmosphere configured to exchange fluxes with a sea ice model and 
and ocean model defined in ClimaOcean
"""
abstract type AbstractAtmosphere end




end