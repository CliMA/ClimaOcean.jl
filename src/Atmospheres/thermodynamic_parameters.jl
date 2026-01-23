import Thermodynamics.Parameters:
    Rv_over_Rd,     # Ratio of the specific gas constants of water vapor over dry air
    R_v,            # Specific gas constant for water vapor
    R_d,            # Specific gas constant for dry air
    kappa_d,        # Ideal gas adiabatic exponent for dry air
    T_0,            # Enthalpy reference temperature
    LH_v0,          # Vaporization enthalpy at the reference temperature
    LH_s0,          # Sublimation enthalpy at the reference temperature
    LH_f0,          # Fusion enthalpy at the reference temperature
    cp_d,           # Heat capacity of dry air at constant pressure
    cp_v,           # Isobaric specific heat capacity of gaseous water vapor
    cp_l,           # Isobaric specific heat capacity of liquid water
    cp_i,           # Isobaric specific heat capacity of water ice
    cv_v,           # Heat capacity of dry air at constant volume
    cv_l,           # Isobaric specific heat capacity of liquid water
    cv_i,           # Isobaric specific heat capacity of liquid water
    e_int_v0,       # what? something about reference internal energy of water vapor
    T_freeze,       # Freezing temperature of _pure_ water
    T_triple,       # Triple point temperature of _pure_ water
    press_triple,   # Triple point pressure of pure water
    T_icenuc,       # Lower temperature limit for the presence of liquid condensate
                    # (below which homogeneous ice nucleation occurs)
    pow_icenuc      # "Power parameter" that controls liquid/ice condensate partitioning
                    # during partial ice nucleation

import ..OceanSeaIceModels: update_net_fluxes!

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
    ConstitutiveParameters(FT = Oceananigans.defaults.FloatType;
                           gas_constant       = 8.3144598,
                           dry_air_molar_mass = 0.02897,
                           water_molar_mass   = 0.018015)

Construct a set of parameters that define the density of moist air,

```math
ρ = p / Rᵐ(q) T,
```

where ``p`` is pressure, ``T`` is temperature, ``q`` defines the partition
of total mass into vapor, liquid, and ice mass fractions, and
``Rᵐ`` is the effective specific gas constant for the mixture,

```math
Rᵐ(q) =
```

where

For more information see [reference docs].
"""
function ConstitutiveParameters(FT = Oceananigans.defaults.FloatType;
                                gas_constant       = 8.3144598,
                                dry_air_molar_mass = 0.02897,
                                water_molar_mass   = 0.018015)

    return ConstitutiveParameters{FT}(convert(FT, gas_constant),
                                      convert(FT, dry_air_molar_mass),
                                      convert(FT, water_molar_mass))
end

const CP{FT} = ConstitutiveParameters{FT} where FT

@inline gas_constant(p::CP)   = p.gas_constant
@inline molmass_dryair(p::CP) = p.dry_air_molar_mass
@inline molmass_water(p::CP)  = p.water_molar_mass
@inline Rv_over_Rd(p::CP)     = molmass_dryair(p) / molmass_water(p)
@inline R_v(p::CP)            = gas_constant(p) / molmass_water(p)
@inline R_d(p::CP)            = gas_constant(p) / molmass_dryair(p)

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
    HeatCapacityParameters(FT = Oceananigans.defaults.FloatType;
                           dry_air_adiabatic_exponent = 2/7,
                           water_vapor_heat_capacity = 1859,
                           liquid_water_heat_capacity = 4181,
                           water_ice_heat_capacity = 2100)

Isobaric heat capacities.
"""
function HeatCapacityParameters(FT = Oceananigans.defaults.FloatType;
                                dry_air_adiabatic_exponent = 2/7,
                                water_vapor_heat_capacity = 1859,
                                liquid_water_heat_capacity = 4181,
                                water_ice_heat_capacity = 2100)

    return HeatCapacityParameters{FT}(convert(FT, dry_air_adiabatic_exponent),
                                      convert(FT, water_vapor_heat_capacity),
                                      convert(FT, liquid_water_heat_capacity),
                                      convert(FT, water_ice_heat_capacity))
end

const HCP{FT} = HeatCapacityParameters{FT} where FT
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

function PhaseTransitionParameters(FT = Oceananigans.defaults.FloatType;
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

const PTP{FT} = PhaseTransitionParameters{FT} where FT
@inline LH_v0(p::PTP)        = p.reference_vaporization_enthalpy
@inline LH_s0(p::PTP)        = p.reference_sublimation_enthalpy
@inline LH_f0(p::PTP)        = LH_s0(p) - LH_v0(p)
@inline T_freeze(p::PTP)     = p.water_freezing_temperature
@inline T_triple(p::PTP)     = p.triple_point_temperature
@inline T_icenuc(p::PTP)     = p.total_ice_nucleation_temperature
@inline pow_icenuc(p::PTP)   = convert(eltype(p), 1) # we shouldn't have the need to set this
@inline press_triple(p::PTP) = p.triple_point_pressure
@inline T_0(p::PTP)          = p.reference_temperature

struct AtmosphereThermodynamicsParameters{FT} <: AbstractThermodynamicsParameters{FT}
    constitutive      :: ConstitutiveParameters{FT}
    heat_capacity     :: HeatCapacityParameters{FT}
    phase_transitions :: PhaseTransitionParameters{FT}
end

const ATP{FT} = AtmosphereThermodynamicsParameters{FT} where FT

Base.eltype(::ATP{FT}) where FT = FT
Base.eltype(::CP{FT})  where FT = FT
Base.eltype(::HCP{FT}) where FT = FT
Base.eltype(::PTP{FT}) where FT = FT

Base.summary(::ATP{FT}) where FT = "AtmosphereThermodynamicsParameters{$FT}"

function Base.show(io::IO, p::AtmosphereThermodynamicsParameters)
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

function AtmosphereThermodynamicsParameters(FT = Oceananigans.defaults.FloatType;
                                            constitutive = ConstitutiveParameters(FT),
                                            phase_transitions = PhaseTransitionParameters(FT),
                                            heat_capacity = HeatCapacityParameters(FT))

    return AtmosphereThermodynamicsParameters(constitutive, heat_capacity, phase_transitions)
end

const ATP = AtmosphereThermodynamicsParameters

@inline R_d(p::ATP)            = R_d(p.constitutive)
@inline R_v(p::ATP)            = R_v(p.constitutive)
@inline gas_constant(p::ATP)   = gas_constant(p.constitutive)
@inline molmass_dryair(p::ATP) = molmass_dryair(p.constitutive)
@inline molmass_water(p::ATP)  = molmass_water(p.constitutive)
@inline Rv_over_Rd(p::ATP)     = Rv_over_Rd(p.constitutive)
@inline LH_v0(p::ATP)          = LH_v0(p.phase_transitions)
@inline LH_s0(p::ATP)          = LH_s0(p.phase_transitions)
@inline LH_f0(p::ATP)          = LH_f0(p.phase_transitions)
@inline T_freeze(p::ATP)       = T_freeze(p.phase_transitions)
@inline T_triple(p::ATP)       = T_triple(p.phase_transitions)
@inline T_icenuc(p::ATP)       = T_icenuc(p.phase_transitions)
@inline pow_icenuc(p::ATP)     = pow_icenuc(p.phase_transitions)
@inline press_triple(p::ATP)   = press_triple(p.phase_transitions)
@inline T_0(p::ATP)            = T_0(p.phase_transitions)

@inline e_int_v0(p::ATP)       = LH_v0(p) - R_v(p) * T_0(p)

@inline cp_v(p::ATP)           = cp_v(p.heat_capacity)
@inline cp_l(p::ATP)           = cp_l(p.heat_capacity)
@inline cp_i(p::ATP)           = cp_i(p.heat_capacity)

@inline cv_l(p::ATP)           = cv_l(p.heat_capacity)
@inline cv_i(p::ATP)           = cv_i(p.heat_capacity)

@inline kappa_d(p::ATP)        = kappa_d(p.heat_capacity)
@inline cp_d(p::ATP)           = R_d(p) / kappa_d(p)
@inline cv_d(p::ATP)           = cp_d(p) - R_d(p)
@inline cv_v(p::ATP)           = cp_v(p) - R_v(p)
