struct PrescribedAtmosphereThermodynamicsParameters{FT} <: AbstractThermodynamicsParameters{FT}
    constitutive      :: ConstitutiveParameters{FT}
    heat_capacity     :: HeatCapacityParameters{FT}
    phase_transitions :: PhaseTransitionParameters{FT}
end

const PATP{FT} = PrescribedAtmosphereThermodynamicsParameters{FT} where FT

Base.eltype(::PATP{FT}) where FT = FT
Base.eltype(::CP{FT})   where FT = FT
Base.eltype(::HCP{FT})  where FT = FT
Base.eltype(::PTP{FT})  where FT = FT

Base.summary(::PATP{FT}) where FT = "PrescribedAtmosphereThermodynamicsParameters{$FT}"

function Base.show(io::IO, p::PrescribedAtmosphereThermodynamicsParameters)
    FT = eltype(p)

    cp = p.constitutive 
    hc = p.heat_capacity
    pt = p.phase_transitions

    return print(io, summary(p), ':', '\n',
        "├── ConstitutiveParameters{$FT}:", '\n',
        "│   ├── gas_constant (R):                      ", prettysummary(cp.gas_constant), '\n',
        "│   ├── dry_air_molar_mass (Mᵈ):               ", prettysummary(cp.dry_air_molar_mass), '\n',
        "│   └── water_molar_mass (Mᵛ):                 ", prettysummary(cp.water_molar_mass), '\n',
        "├── HeatCapacityParameters{$FT}:", '\n',
        "│   ├── dry_air_adiabatic_exponent (κᵈ):       ", prettysummary(hc.dry_air_adiabatic_exponent), '\n',
        "│   ├── water_vapor_heat_capacity (cᵖᵛ):       ", prettysummary(hc.water_vapor_heat_capacity), '\n',
        "│   ├── liquid_water_heat_capacity (cᵖˡ):      ", prettysummary(hc.liquid_water_heat_capacity), '\n',
        "│   └── water_ice_heat_capacity (cᵖⁱ):         ", prettysummary(hc.water_ice_heat_capacity), '\n',
        "└── PhaseTransitionParameters{$FT}", '\n',
        "    ├── reference_vaporization_enthalpy (ℒᵛ⁰): ", prettysummary(pt.reference_vaporization_enthalpy), '\n',
        "    ├── reference_sublimation_enthalpy  (ℒˢ⁰): ", prettysummary(pt.reference_sublimation_enthalpy), '\n',
        "    ├── reference_temperature (T⁰):            ", prettysummary(pt.reference_temperature), '\n',    
        "    ├── triple_point_temperature (Tᵗʳ):        ", prettysummary(pt.triple_point_temperature), '\n',
        "    ├── triple_point_pressure (pᵗʳ):           ", prettysummary(pt.triple_point_pressure), '\n',   
        "    ├── water_freezing_temperature (Tᶠ):       ", prettysummary(pt.water_freezing_temperature), '\n',
        "    └── total_ice_nucleation_temperature (Tⁱ): ", prettysummary(pt.total_ice_nucleation_temperature))
end

function PrescribedAtmosphereThermodynamicsParameters(FT = Float64;
                                                      constitutive = ConstitutiveParameters(FT),
                                                      phase_transitions = PhaseTransitionParameters(FT),
                                                      heat_capacity = HeatCapacityParameters(FT))

    return PrescribedAtmosphereThermodynamicsParameters(constitutive, heat_capacity, phase_transitions)
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

@inline R_d(p::PATP)            = R_d(p.constitutive)
@inline R_v(p::PATP)            = R_v(p.constitutive)
@inline gas_constant(p::PATP)   = gas_constant(p.constitutive)
@inline molmass_dryair(p::PATP) = molmass_dryair(p.constitutive)
@inline molmass_water(p::PATP)  = molmass_water(p.constitutive)
@inline molmass_ratio(p::PATP)  = molmass_ratio(p.constitutive)
@inline LH_v0(p::PATP)          = LH_v0(p.phase_transitions)
@inline LH_s0(p::PATP)          = LH_s0(p.phase_transitions)
@inline LH_f0(p::PATP)          = LH_f0(p.phase_transitions)
@inline T_freeze(p::PATP)       = T_freeze(p.phase_transitions)
@inline T_triple(p::PATP)       = T_triple(p.phase_transitions)
@inline T_icenuc(p::PATP)       = T_icenuc(p.phase_transitions)
@inline pow_icenuc(p::PATP)     = pow_icenuc(p.phase_transitions)
@inline press_triple(p::PATP)   = press_triple(p.phase_transitions)
@inline T_0(p::PATP)            = T_0(p.phase_transitions)

@inline e_int_v0(p::PATP)       = LH_v0(p) - R_v(p) * T_0(p)

@inline cp_v(p::PATP)           = cp_v(p.heat_capacity)
@inline cp_l(p::PATP)           = cp_l(p.heat_capacity)
@inline cp_i(p::PATP)           = cp_i(p.heat_capacity)

@inline cv_l(p::PATP)           = cv_l(p.heat_capacity)
@inline cv_i(p::PATP)           = cv_i(p.heat_capacity)

@inline kappa_d(p::PATP)        = kappa_d(p.heat_capacity)
@inline cp_d(p::PATP)           = R_d(p) / kappa_d(p)
@inline cv_d(p::PATP)           = cp_d(p) - R_d(p)
@inline cv_v(p::PATP)           = cp_v(p) - R_v(p)

#####
##### Prescribed atmosphere (as opposed to dynamically evolving / prognostic)
#####

struct PrescribedAtmosphere{FT, M, G, U, P, C, F, I, R, TP, TI} 
    grid :: G
    metadata :: M
    velocities :: U
    pressure :: P
    tracers :: C
    freshwater_flux :: F
    auxiliary_freshwater_flux :: I
    downwelling_radiation :: R
    thermodynamics_parameters :: TP
    times :: TI
    reference_height :: FT
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
    print(io, "├── reference_height: ", prettysummary(pa.reference_height), '\n')
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

function default_atmosphere_pressure(grid, times)
    pa = FieldTimeSeries{Center, Center, Nothing}(grid, times)
    parent(pa) .= 101325
    return pa
end

"""
    PrescribedAtmosphere(grid, times;
                         metadata = nothing,
                         reference_height = 10, # meters
                         boundary_layer_height = 600 # meters,
                         thermodynamics_parameters = PrescribedAtmosphereThermodynamicsParameters(FT),
                         auxiliary_freshwater_flux = nothing,
                         velocities            = default_atmosphere_velocities(grid, times),
                         tracers               = default_atmosphere_tracers(grid, times),
                         pressure              = default_atmosphere_pressure(grid, times),
                         freshwater_flux       = default_freshwater_flux(grid, times),
                         downwelling_radiation = default_downwelling_radiation(grid, times))

Return a representation of a prescribed time-evolving atmospheric
state with data given at `times`.
"""
function PrescribedAtmosphere(grid, times;
                              metadata = nothing,  
                              reference_height = convert(eltype(grid), 10),
                              boundary_layer_height = convert(eltype(grid), 600),
                              thermodynamics_parameters = nothing,
                              auxiliary_freshwater_flux = nothing,
                              velocities            = default_atmosphere_velocities(grid, times),
                              tracers               = default_atmosphere_tracers(grid, times),
                              pressure              = default_atmosphere_pressure(grid, times),
                              freshwater_flux       = default_freshwater_flux(grid, times),
                              downwelling_radiation = default_downwelling_radiation(grid, times))

    FT = eltype(grid)
    if isnothing(thermodynamics_parameters)
        thermodynamics_parameters = PrescribedAtmosphereThermodynamicsParameters(FT)
    end

    return PrescribedAtmosphere(grid,
                                metadata,
                                velocities,
                                pressure,
                                tracers,
                                freshwater_flux,
                                auxiliary_freshwater_flux,
                                downwelling_radiation,
                                thermodynamics_parameters,
                                times,
                                convert(FT, reference_height),
                                convert(FT, boundary_layer_height))
end

# No need to timestep!
time_step!(atmos::PrescribedAtmosphere) = nothing

update_model_field_time_series!(::Nothing, time) = nothing

function update_model_field_time_series!(atmos::PrescribedAtmosphere, time)
    ftses = extract_field_time_series(atmos)
    for fts in ftses
        update_field_time_series!(fts, time)
    end

    return nothing
end

struct TwoBandDownwellingRadiation{SW, LW}
    shortwave :: SW
    longwave :: LW
end

"""
    TwoBandDownwellingRadiation(shortwave=nothing, longwave=nothing)

Return a two-band model for downwelling radiation (split in a shortwave band
and a longwave band) that passes through the atmosphere and arrives at the surface of ocean
or sea ice.
"""
TwoBandDownwellingRadiation(; shortwave=nothing, longwave=nothing) =
    TwoBandDownwellingRadiation(shortwave, longwave)

Adapt.adapt_structure(to, tsdr::TwoBandDownwellingRadiation) =
    TwoBandDownwellingRadiation(adapt(to, tsdr.shortwave),
                                adapt(to, tsdr.longwave))

# A prescribed atmosphere does not need fluxes!
regrid_fluxes_to_atmospheric_model!(atmos::PrescribedAtmosphere, args...) = nothing

function interpolate_atmospheric_state!(surface_atmosphere_state,
                                        interpolated_prescribed_freshwater_flux,
                                        atmosphere::PrescribedAtmosphere, 
                                        ocean_grid, clock)

    atmosphere_grid = atmosphere.grid

    # We use .data here to save parameter space (unlike Field, adapt_structure for
    # fts = FieldTimeSeries does not return fts.data)
    atmosphere_velocities = (u = atmosphere.velocities.u.data,
                             v = atmosphere.velocities.v.data)

    atmosphere_tracers = (T = atmosphere.tracers.T.data,
                          q = atmosphere.tracers.q.data)

    Qs = atmosphere.downwelling_radiation.shortwave
    Qℓ = atmosphere.downwelling_radiation.longwave
    downwelling_radiation = (shortwave=Qs.data, longwave=Qℓ.data)
    freshwater_flux = map(ϕ -> ϕ.data, atmosphere.freshwater_flux)
    atmosphere_pressure = atmosphere.pressure.data

    # Extract info for time-interpolation
    u = atmosphere.velocities.u # for example
    atmosphere_times = u.times
    atmosphere_backend = u.backend
    atmosphere_time_indexing = u.time_indexing

    Nx, Ny, Nz = size(ocean_grid)
    single_column_grid = Nx == 1 && Ny == 1

    if single_column_grid
        kernel_parameters = KernelParameters(1:1, 1:1)
    else
        # Compute fluxes into halo regions, ie from 0:Nx+1 and 0:Ny+1
        kernel_parameters = KernelParameters(0:Nx+1, 0:Ny+1)
    end    
    
    launch!(architecture(ocean_grid), ocean_grid, kernel_parameters,
            _interpolate_primary_atmospheric_state!,
            surface_atmosphere_state,
            ocean_grid,
            clock,
            atmosphere_velocities,
            atmosphere_tracers,
            atmosphere_pressure,
            downwelling_radiation,
            freshwater_flux,
            atmosphere_grid,
            atmosphere_times,
            atmosphere_backend,
            atmosphere_time_indexing)

    # Separately interpolate the auxiliary freshwater fluxes, which may
    # live on a different grid than the primary fluxes and atmospheric state.
    
    auxiliary_freshwater_flux = atmosphere.auxiliary_freshwater_flux

    if !isnothing(auxiliary_freshwater_flux)
        # TODO: do not assume that `auxiliary_freshater_flux` is a tuple
        auxiliary_data = map(ϕ -> ϕ.data, auxiliary_freshwater_flux)

        first_auxiliary_flux    = first(auxiliary_freshwater_flux)
        auxiliary_grid          = first_auxiliary_flux.grid
        auxiliary_times         = first_auxiliary_flux.times
        auxiliary_backend       = first_auxiliary_flux.backend
        auxiliary_time_indexing = first_auxiliary_flux.time_indexing

        launch!(architecture(ocean_grid), ocean_grid, kernel_parameters,
                _interpolate_auxiliary_freshwater_flux!,
                interpolated_prescribed_freshwater_flux,
                ocean_grid,
                clock,
                auxiliary_data,
                auxiliary_grid,
                auxiliary_times,
                auxiliary_backend,
                auxiliary_time_indexing)
    end

    return nothing
end

@kernel function _interpolate_primary_atmospheric_state!(surface_atmos_state,
                                                         grid,
                                                         clock,
                                                         atmos_velocities,
                                                         atmos_tracers,
                                                         atmos_pressure,
                                                         downwelling_radiation,
                                                         prescribed_freshwater_flux,
                                                         atmos_grid,
                                                         atmos_times,
                                                         atmos_backend,
                                                         atmos_time_indexing)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) # index of the top ocean cell

    @inbounds begin
        # Atmos state, which is _assumed_ to exist at location = (c, c, nothing)
        # The third index "k" should not matter but we put the correct index to get
        # a surface node anyways.
        atmos_args = (atmos_grid, atmos_times, atmos_backend, atmos_time_indexing)
        X = _node(i, j, kᴺ + 1, grid, c, c, f)
        time = Time(clock.time)

        uₐ = interp_atmos_time_series(atmos_velocities.u, X, time, atmos_args...)
        vₐ = interp_atmos_time_series(atmos_velocities.v, X, time, atmos_args...)
        Tₐ = interp_atmos_time_series(atmos_tracers.T,    X, time, atmos_args...)
        qₐ = interp_atmos_time_series(atmos_tracers.q,    X, time, atmos_args...)
        pₐ = interp_atmos_time_series(atmos_pressure,     X, time, atmos_args...)

        Qs = interp_atmos_time_series(downwelling_radiation.shortwave, X, time, atmos_args...)
        Qℓ = interp_atmos_time_series(downwelling_radiation.longwave,  X, time, atmos_args...)

        # Usually precipitation
        Mh = interp_atmos_time_series(prescribed_freshwater_flux, X, time, atmos_args...)

        # Convert atmosphere velocities (defined on a latitude-longitude grid) to 
        # the frame of reference of the native grid
        uₐ, vₐ = intrinsic_vector(i, j, kᴺ, grid, uₐ, vₐ)
    
        surface_atmos_state.u[i, j, 1] = uₐ
        surface_atmos_state.v[i, j, 1] = vₐ
        surface_atmos_state.T[i, j, 1] = Tₐ
        surface_atmos_state.p[i, j, 1] = pₐ
        surface_atmos_state.q[i, j, 1] = qₐ
        surface_atmos_state.Qs[i, j, 1] = Qs
        surface_atmos_state.Qℓ[i, j, 1] = Qℓ
        surface_atmos_state.Mp[i, j, 1] = Mh
    end
end

@kernel function _interpolate_auxiliary_freshwater_flux!(freshwater_flux,
                                                         grid,
                                                         clock,
                                                         auxiliary_freshwater_flux,
                                                         auxiliary_grid,
                                                         auxiliary_times,
                                                         auxiliary_backend,
                                                         auxiliary_time_indexing)

    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3) # index of the top ocean cell

    @inbounds begin
        X = _node(i, j, kᴺ + 1, grid, c, c, f)
        time = Time(clock.time)
        Mr = interp_atmos_time_series(auxiliary_freshwater_flux, X, time,
                                      auxiliary_grid,
                                      auxiliary_times,
                                      auxiliary_backend,
                                      auxiliary_time_indexing)

        freshwater_flux[i, j, 1] += Mr
    end
end
