module PrescribedAtmospheres

using Oceananigans.Utils: prettysummary

using Adapt
using Thermodynamics.Parameters: AbstractThermodynamicsParameters

import Thermodynamics.Parameters:
    gas_constant,   #
    molmass_dryair, # Molar mass of dry air (without moisture)
    molmass_water,  # Molar mass of gaseous water vapor
    molmass_ratio,  # Ratio of the molar masses of dry air to water vapor
    R_v,            # Specific gas constant for water vapor
    R_d,            # Specific gas constant for dry air
    kappa_d,        # Ideal gas adiabatic exponent for dry air
    T_0,            # Enthalpy reference temperature
    LH_v0,          # Vaporization enthalpy at the reference temperature
    LH_s0,          # Sublimation enthalpy at the reference temperature
    LH_f0,          # Fusionn enthalpy at the reference temperature
    cp_d,           # Heat capacity of dry air at constant pressure
    cp_v,           # Isobaric specific heat capacity of gaseous water vapor
    cp_l,           # Isobaric specific heat capacity of liquid water
    cp_i,           # Isobaric specific heat capacity of water ice
    cv_v,           # Heat capacity of dry air at constant volume
    cv_l,           # Isobaric specific heat capacity of liquid water
    cv_i,           # Isobaric specific heat capacity of liquid water
    e_int_v0,       # what? someting about reference internal energy of water vapor
    T_freeze,       # Freezing temperature of _pure_ water
    T_triple,       # Triple point temperature of _pure_ water
    press_triple,   # Triple point pressure of pure water
    T_icenuc,       # Lower temperature limit for the presence of liquid condensate
                    # (below which homogeneous ice nucleation occurs)
    pow_icenuc      # "Power parameter" that controls liquid/ice condensate partitioning
                    # during partial ice nucleation

import ..OceanSeaIceModels:
    surface_velocities,
    surface_tracers,
    downwelling_radiation,
    freshwater_flux

#####
##### Atmospheric thermodynamics parameters
#####

struct ConstitutiveParameters{FT} <: AbstractThermodynamicsParameters{FT}
    gas_constant       :: FT
    dry_air_molar_mass :: FT
    water_molar_mass   :: FT
end

function Base.summary(p::ConstitutiveParameters{FT}) where FT
    return string("ConstitutiveParameters{$FT}(",
                    "R=", prettysummary(p.gas_constant),
                  ", Mᵈ=", prettysummary(p.dry_air_molar_mass),
                  ", Mᵛ=", prettysummary(p.water_molar_mass), ")")
end

Base.show(io::IO, p::ConstitutiveParameters) = print(io, summary(p))

"""
    ConstitutiveParameters(FT; gas_constant       = 8.3144598,
                               dry_air_molar_mass = 0.02897,
                               water_molar_mass   = 0.018015)

Construct a set of parameters that define the density of moist air,

```math
ρ = p / Rᵐ(q) T,
```

where ``p`` is pressure, ``T`` is temperature, ``q`` defines the partition
of total mass into vapor, liqiud, and ice mass fractions, and
``Rᵐ`` is the effective specific gas constant for the mixture,

```math
Rᵐ(q) =
```

where

For more information see [reference docs].
"""
function ConstitutiveParameters(FT = Float64;
                                gas_constant       = 8.3144598,
                                dry_air_molar_mass = 0.02897,
                                water_molar_mass   = 0.018015)

    return ConstitutiveParameters{FT}(convert(FT, gas_constant),
                                      convert(FT, dry_air_molar_mass),
                                      convert(FT, water_molar_mass))
end

const CP = ConstitutiveParameters

@inline gas_constant(p::CP)     = p.gas_constant
@inline molmass_dryair(p::CP)   = p.dry_air_molar_mass
@inline molmass_water(p::CP)    = p.water_molar_mass
@inline molmass_ratio(p::CP)    = molmass_dryair(p) / molmass_water(p)
@inline R_v(p::CP)              = gas_constant(p) / molmass_water(p)
@inline R_d(p::CP)              = gas_constant(p) / molmass_dryair(p)

struct HeatCapacityParameters{FT} <: AbstractThermodynamicsParameters{FT}
    dry_air_adiabatic_exponent :: FT
    water_vapor_heat_capacity  :: FT
    liquid_water_heat_capacity :: FT
    water_ice_heat_capacity    :: FT
end

function Base.summary(p::HeatCapacityParameters{FT}) where FT
    return string("HeatCapacityParameters{$FT}(",
                    "κᵈ=", prettysummary(p.dry_air_adiabatic_exponent),
                  ", cᵖᵛ=", prettysummary(p.water_vapor_heat_capacity),
                  ", cᵖˡ=", prettysummary(p.liquid_water_heat_capacity),
                  ", cᵖⁱ=", prettysummary(p.water_ice_heat_capacity))
end

Base.show(io::IO, p::HeatCapacityParameters) = print(io, summary(p))

"""
    HeatCapacityParameters(FT = Float64,
                           dry_air_adiabatic_exponent = 2/7,
                           water_vapor_heat_capacity = 1859,
                           liquid_water_heat_capacity = 4181,
                           water_ice_heat_capacity = 2100)

Isobaric heat capacities.
"""
function HeatCapacityParameters(FT = Float64;
                                dry_air_adiabatic_exponent = 2/7,
                                water_vapor_heat_capacity = 1859,
                                liquid_water_heat_capacity = 4181,
                                water_ice_heat_capacity = 2100)

    return HeatCapacityParameters{FT}(convert(FT, dry_air_adiabatic_exponent),
                                      convert(FT, water_vapor_heat_capacity),
                                      convert(FT, liquid_water_heat_capacity),
                                      convert(FT, water_ice_heat_capacity))
end

const HCP = HeatCapacityParameters
@inline cp_v(p::HCP)    = p.water_vapor_heat_capacity
@inline cp_l(p::HCP)    = p.liquid_water_heat_capacity
@inline cp_i(p::HCP)    = p.water_ice_heat_capacity
@inline cv_l(p::HCP)    = cp_l(p)
@inline cv_i(p::HCP)    = cp_i(p)
@inline kappa_d(p::HCP) = p.dry_air_adiabatic_exponent

struct PhaseTransitionParameters{FT} <: AbstractThermodynamicsParameters{FT}
    reference_vaporization_enthalpy  :: FT
    reference_sublimation_enthalpy   :: FT
    reference_temperature            :: FT
    triple_point_temperature         :: FT
    triple_point_pressure            :: FT
    water_freezing_temperature       :: FT
    total_ice_nucleation_temperature :: FT
end

function Base.summary(p::PhaseTransitionParameters{FT}) where FT
    return string("PhaseTransitionParameters{$FT}(",
                    "ℒᵛ⁰=", prettysummary(p.reference_vaporization_enthalpy),
                  ", ℒˢ⁰=", prettysummary(p.reference_sublimation_enthalpy),
                  ", T⁰=", prettysummary(p.reference_temperature),
                  ", Tᵗʳ=", prettysummary(p.triple_point_temperature),
                  ", pᵗʳ=", prettysummary(p.triple_point_pressure),
                  ", Tᶠ=", prettysummary(p.water_freezing_temperature),
                  ", Tⁱⁿ=", prettysummary(p.total_ice_nucleation_temperature), ')')
end

Base.show(io::IO, p::PhaseTransitionParameters) = print(io, summary(p))

function PhaseTransitionParameters(FT = Float64;
                                   reference_vaporization_enthalpy = 2500800,
                                   reference_sublimation_enthalpy = 2834400,
                                   reference_temperature = 273.16,
                                   triple_point_temperature = 273.16,
                                   triple_point_pressure = 611.657,
                                   water_freezing_temperature = 273.15,
                                   total_ice_nucleation_temperature = 233)

    return PhaseTransitionParameters{FT}(convert(FT, reference_vaporization_enthalpy),
                                         convert(FT, reference_sublimation_enthalpy),
                                         convert(FT, reference_temperature),
                                         convert(FT, triple_point_temperature),
                                         convert(FT, triple_point_pressure),
                                         convert(FT, water_freezing_temperature),
                                         convert(FT, total_ice_nucleation_temperature))
end

const PTP = PhaseTransitionParameters
@inline LH_v0(p::PTP)        = p.reference_vaporization_enthalpy
@inline LH_s0(p::PTP)        = p.reference_sublimation_enthalpy
@inline LH_f0(p::PTP)        = LH_s0(p) - LH_v0(p)
@inline T_freeze(p::PTP)     = p.water_freezing_temperature
@inline T_triple(p::PTP)     = p.triple_point_temperature
@inline T_icenuc(p::PTP)     = p.total_ice_nucleation_temperature
@inline pow_icenuc(p::PTP)   = convert(eltype(p), 1) # we shouldn't have the need to set this
@inline press_triple(p::PTP) = p.triple_point_pressure
@inline T_0(p::PTP)          = p.reference_temperature

struct PrescribedAtmosphereThermodynamicsParameters{FT} <: AbstractThermodynamicsParameters{FT}
    constitutive      :: ConstitutiveParameters{FT}
    heat_capacity     :: HeatCapacityParameters{FT}
    phase_transitions :: PhaseTransitionParameters{FT}
end

const PATP{FT} = PrescribedAtmosphereThermodynamicsParameters{FT} where FT

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

struct PrescribedAtmosphere{G, U, P, C, F, R, TP, TI, FT}
    grid :: G
    velocities :: U
    pressure :: P
    tracers :: C
    freshwater_flux :: F
    downwelling_radiation :: R
    thermodynamics_parameters :: TP
    times :: TI
    reference_height :: FT
end

Base.summary(::PrescribedAtmosphere) = "PrescribedAtmosphere"
Base.show(io::IO, pa::PrescribedAtmosphere) = print(io, summary(pa))

"""
    PrescribedAtmosphere(times;
                         reference_height,
                         velocities = nothing,
                         pressure = nothing,
                         freshwater_flux = nothing,
                         downwelling_radiation = nothing,
                         tracers = nothing)

Return a representation of a prescribed time-evolving atmospheric
state with data given at `times`.
"""
function PrescribedAtmosphere(times, FT=Float64;
                              reference_height,
                              velocities = nothing,
                              pressure = nothing,
                              freshwater_flux = nothing,
                              downwelling_radiation = nothing,
                              thermodynamics_parameters = PrescribedAtmosphereThermodynamicsParameters(FT),
                              grid = nothing,
                              tracers = nothing)

    if isnothing(grid) # try to find it
        u = first(velocities)
        grid = u.grid
    end

    return PrescribedAtmosphere(grid,
                                velocities,
                                pressure,
                                tracers,
                                freshwater_flux,
                                downwelling_radiation,
                                thermodynamics_parameters,
                                times,
                                convert(FT, reference_height))
end

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

Adapt.adapt_structure(to, tsdr::TwoStreamDownwellingRadiation) =
    TwoStreamDownwellingRadiation(adapt(to, tsdr.shortwave),
                                  adapt(to, tsdr.longwave))

end # module

