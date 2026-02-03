#####
##### Prescribed atmosphere (as opposed to dynamically evolving / prognostic)
#####

mutable struct PrescribedAtmosphere{FT, G, T, U, P, C, F, I, R, TP, TI}
    grid :: G
    clock :: Clock{T}
    velocities :: U
    pressure :: P
    tracers :: C
    freshwater_flux :: F
    auxiliary_freshwater_flux :: I
    downwelling_radiation :: R
    thermodynamics_parameters :: TP
    times :: TI
    surface_layer_height :: FT
    boundary_layer_height :: FT
end

function Base.summary(pa::PrescribedAtmosphere{FT}) where FT
    Nx, Ny, Nz = size(pa.grid)
    Nt = length(pa.times)
    sz_str = string(Nx, "×", Ny, "×", Nz, "×", Nt)
    return string(sz_str, " PrescribedAtmosphere{$FT}")
end

function Base.show(io::IO, pa::PrescribedAtmosphere)
    print(io, summary(pa), " on ", grid_name(pa.grid), ":", '\n')
    print(io, "├── times: ", prettysummary(pa.times), '\n')
    print(io, "├── surface_layer_height: ", prettysummary(pa.surface_layer_height), '\n')
    print(io, "└── boundary_layer_height: ", prettysummary(pa.boundary_layer_height))
end

function default_atmosphere_velocities(grid, times)
    ua = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    va = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    return (u=ua, v=va)
end

function default_atmosphere_tracers(grid, times)
    Ta = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    qa = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    parent(Ta) .= 273.15 + 20
    return (T=Ta, q=qa)
end

function default_downwelling_radiation(grid, times)
    Qℓ = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    Qs = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    return TwoBandDownwellingRadiation(shortwave=Qs, longwave=Qℓ)
end

function default_freshwater_flux(grid, times)
    rain = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    snow = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    return (; rain, snow)
end

""" The standard unit of atmospheric pressure; 1 standard atmosphere (atm) = 101,325 Pascals (Pa)
in SI units. This is approximately equal to the mean sea-level atmospheric pressure on Earth. """
function default_atmosphere_pressure(grid, times)
    pa = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    parent(pa) .= 101325
    return pa
end

@inline function update_state!(atmos::PrescribedAtmosphere)
    time = Time(atmos.clock.time)
    ftses = extract_field_time_series(atmos)

    for fts in ftses
        update_field_time_series!(fts, time)
    end
    return nothing
end

@inline function time_step!(atmos::PrescribedAtmosphere, Δt)
    tick!(atmos.clock, Δt)

    update_state!(atmos)

    return nothing
end

@inline thermodynamics_parameters(atmos::Nothing) = nothing
@inline thermodynamics_parameters(atmos::PrescribedAtmosphere) = atmos.thermodynamics_parameters
@inline surface_layer_height(atmos::PrescribedAtmosphere) = atmos.surface_layer_height
@inline boundary_layer_height(atmos::PrescribedAtmosphere) = atmos.boundary_layer_height

# No need to compute anything here...
update_net_fluxes!(coupled_model, ::PrescribedAtmosphere) = nothing

"""
    PrescribedAtmosphere(grid, times=[zero(grid)];
                         clock = Clock{Float64}(time = 0),
                         surface_layer_height = 10, # meters
                         boundary_layer_height = 512 # meters,
                         thermodynamics_parameters = AtmosphereThermodynamicsParameters(eltype(grid)),
                         auxiliary_freshwater_flux = nothing,
                         velocities            = default_atmosphere_velocities(grid, times),
                         tracers               = default_atmosphere_tracers(grid, times),
                         pressure              = default_atmosphere_pressure(grid, times),
                         freshwater_flux       = default_freshwater_flux(grid, times),
                         downwelling_radiation = default_downwelling_radiation(grid, times))

Return a representation of a prescribed time-evolving atmospheric
state with data given at `times`.
"""
function PrescribedAtmosphere(grid, times=[zero(grid)];
                              clock = Clock{Float64}(time = 0),
                              surface_layer_height = 10,
                              boundary_layer_height = 512,
                              thermodynamics_parameters = AtmosphereThermodynamicsParameters(eltype(grid)),
                              auxiliary_freshwater_flux = nothing,
                              velocities            = default_atmosphere_velocities(grid, times),
                              tracers               = default_atmosphere_tracers(grid, times),
                              pressure              = default_atmosphere_pressure(grid, times),
                              freshwater_flux       = default_freshwater_flux(grid, times),
                              downwelling_radiation = default_downwelling_radiation(grid, times))

    FT = eltype(grid)
    if isnothing(thermodynamics_parameters)
        thermodynamics_parameters = AtmosphereThermodynamicsParameters(FT)
    end

    atmosphere = PrescribedAtmosphere(grid,
                                      clock,
                                      velocities,
                                      pressure,
                                      tracers,
                                      freshwater_flux,
                                      auxiliary_freshwater_flux,
                                      downwelling_radiation,
                                      thermodynamics_parameters,
                                      times,
                                      convert(FT, surface_layer_height),
                                      convert(FT, boundary_layer_height))
    update_state!(atmosphere)

    return atmosphere
end

struct TwoBandDownwellingRadiation{SW, LW}
    shortwave :: SW
    longwave :: LW
end

"""
    TwoBandDownwellingRadiation(shortwave=nothing, longwave=nothing)

Return a two-band model for downwelling radiation (split into a shortwave band
and a longwave band) that passes through the atmosphere and arrives at the surface of ocean
or sea ice.
"""
TwoBandDownwellingRadiation(; shortwave=nothing, longwave=nothing) =
    TwoBandDownwellingRadiation(shortwave, longwave)

Adapt.adapt_structure(to, tsdr::TwoBandDownwellingRadiation) =
    TwoBandDownwellingRadiation(adapt(to, tsdr.shortwave),
                                adapt(to, tsdr.longwave))

#####
##### Chekpointing
#####

import Oceananigans: prognostic_state, restore_prognostic_state!

function prognostic_state(atmos::PrescribedAtmosphere) 
    return (; clock = prognostic_state(atmos.clock))
end

function restore_prognostic_state!(atmos::PrescribedAtmosphere, state) 
    restore_prognostic_state!(atmos.clock, state.clock)
    update_state!(atmos)
    return atmos
end

restore_prognostic_state!(atmos::PrescribedAtmosphere, ::Nothing) = atmos
