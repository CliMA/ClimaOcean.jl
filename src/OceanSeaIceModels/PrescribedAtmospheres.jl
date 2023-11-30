module PrescribedAtmospheres

import ..OceanSeaIceModels: surface_velocities, surface_tracers

struct PrescribedAtmosphere{U, F, R, C, T}
    velocities :: U
    freshwater_flux :: F
    downwelling_radiation :: R
    tracers :: C
    times :: T
end

"""
    PrescribedAtmosphere(times;
                         velocities = nothing,
                         freshwater_flux = nothing,
                         downwelling_radiation = nothing,
                         tracers = nothing)

Return a representation of a prescribed time-evolving atmospheric
state with data given at `times`.
"""
function PrescribedAtmosphere(times;
                              velocities = nothing,
                              freshwater_flux = nothing,
                              downwelling_radiation = nothing,
                              tracers = nothing)

    return PrescribedAtmosphere(velocities,
                                freshwater_flux,
                                downwelling_radiation,
                                tracers,
                                times)
end

surface_velocities(pa::PrescribedAtmosphere) = pa.velocities
surface_tracers(pa::PrescribedAtmosphere) = pa.tracers
freshwater_fluxes(pa::PrescribedAtmosphere) = pa.freshwater_fluxes

struct TwoStreamDownwellingRadiation{SW, LW}
    shortwave :: SW
    longwave :: LW
end

"""
    TwoStreamDownwellingRadiation(shortwave=nothing, longwave=nothing)

Return a two-stream model for downwelling radiation that
passes through he atmosphere and arrives at the surface of ocean
or sea ice.
"""
TwoStreamDownwellingRadiation(; shortwave=nothing, longwave=nothing) =
    TwoStreamDownwellingRadiation(shortwave, longwave)
 
end # module

