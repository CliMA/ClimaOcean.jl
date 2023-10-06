module AtmosphericForcings

# We generally have 2 types of atmospheric forcing: Prescribed fluxes and
# Prescribed atmospheric state (to treat with bulk formulae)

export PrescribedAtmosphere, PrescribedFluxes

abstract type AbstractAtmospericForcing end

struct PrescribedAtmosphere{} <: AbstractAtmospericForcing


end




end