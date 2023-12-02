module PrescribedAtmospheres

import ..OceanSeaIceModels: surface_velocities, surface_tracers, downwelling_radiation, freshwater_flux
import ..OceanSeaIceModels: density, specific_heat

struct PrescribedAtmosphere{U, F, R, C, ρ, CP, T}
    velocities :: U
    freshwater_flux :: F
    downwelling_radiation :: R
    tracers :: C
    density :: ρ
    specific_heat :: CP
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
                              density = 1.2,
                              specific_heat = 1.0005,
                              tracers = nothing)

    return PrescribedAtmosphere(velocities,
                                freshwater_flux,
                                downwelling_radiation,
                                tracers,
                                density,
                                specific_heat,
                                times)
end

surface_velocities(pa::PrescribedAtmosphere) = pa.velocities
surface_tracers(pa::PrescribedAtmosphere) = pa.tracers
freshwater_flux(pa::PrescribedAtmosphere) = pa.freshwater_flux
downwelling_radiation(pa::PrescribedAtmosphere) = pa.downwelling_radiation
density(pa::PrescribedAtmosphere) = pa.density
specific_heat(pa::PrescribedAtmosphere) = pa.specific_heat

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

