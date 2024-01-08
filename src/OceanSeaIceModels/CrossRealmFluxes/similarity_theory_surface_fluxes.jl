using Oceananigans.Utils: prettysummary

using CLIMAParameters
using CLIMAParameters: AliasParamDict

using SurfaceFluxes.Parameters: SurfaceFluxesParameters, AbstractSurfaceFluxesParameters
using SurfaceFluxes.UniversalFunctions: BusingerParams

using Thermodynamics.Parameters: AbstractThermodynamicsParameters

# TODO: write parameter meaning here
import Thermodynamics.Parameters:
    gas_constant,
    molmass_dryair,
    molmass_water,
    kappa_d,
    LH_v0,
    LH_s0,
    cp_v,
    cp_l,
    cp_i,
    T_freeze,
    T_triple,
    T_icenuc,
    press_triple,
    T_0

#####
##### Similarity Theory bulk turbulent fluxes
#####

struct SimilarityTheoryTurbulentFluxes{FT, ΔU, UF, TP} <: AbstractSurfaceFluxesParameters
    gravitational_acceleration :: FT
    von_karman_constant :: FT
    bulk_velocity_scale :: ΔU
    universal_function :: UF
    thermodynamics_parameters :: TP
end

Base.summary(::SimilarityTheoryTurbulentFluxes{FT}) where FT = "SimilarityTheoryTurbulentFluxes{$FT}"

function Base.show(io::IO, fluxes::SimilarityTheoryTurbulentFluxes)
    print(io, summary(fluxes), '\n',
          "├── gravitational_acceleration: ", prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",  prettysummary(fluxes.von_karman_constant), '\n',
          "├── bulk_velocity_scale: ",  summary(fluxes.bulk_velocity_scale), '\n',
          "├── universal_function: ", summary(fluxes.universal_function), '\n',
          "└── thermodynamics_parameters: ", summary(fluxes.thermodynamics_parameters))
end

function SimilarityTheoryTurbulentFluxes(FT = Float64;
                                         gravitational_acceleration = 9.80665,
                                         bulk_velocity_scale = nothing,
                                         von_karman_constant = 0.4,
                                         universal_function = default_universal_function_parameters(FT),
                                         thermodynamics_parameters = HierarchicalThermodynamicsParameters(FT))

    return SimilarityTheoryTurbulentFluxes(gravitational_acceleration,
                                           von_karman_constant,
                                           bulk_velocity_scale,
                                           universal_function,
                                           thermodynamics_parameters)
end

#####
##### Thermodynamics parameters
#####

struct ConstitutiveParameters{FT} <: AbstractThermodynamicsParameters{FT}
    gas_constant               :: FT
    dry_air_molar_mass         :: FT
    water_molar_mass           :: FT
end

"""
    ConstitutiveParameters(FT; gas_constant               = 8.3144598,
                               dry_air_molar_mass         = 0.02897,
                               water_molar_mass           = 0.018015)

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

    return ConstitutiveParameters(convert(FT, gas_constant),
                                  convert(FT, dry_air_molar_mass),
                                  convert(FT, water_molar_mass))
end

const CP = ConstitutiveParameters

gas_constant(p::CP)   = p.gas_constant
molmass_dryair(p::CP) = p.dry_air_molar_mass
molmass_water(p::CP)  = p.water_molar_mass

struct HeatCapacityParameters{FT} <: AbstractThermodynamicsParameters{FT}
    dry_air_adiabatic_exponent :: FT
    water_vapor_heat_capacity  :: FT
    liquid_water_heat_capacity :: FT
    ice_heat_capacity          :: FT
end

function HeatCapacityParameters(FT = Float64;
                                dry_air_adiabatic_exponent = 2/7,
                                water_vapor_heat_capacity = 1859,
                                liquid_water_heat_capacity = 4181,
                                ice_heat_capacity = 2100)

    return HeatCapacityParameters(convert(FT, dry_air_adiabatic_exponent),
                                  convert(FT, water_vapor_heat_capacity),
                                  convert(FT, liquid_water_heat_capacity),
                                  convert(FT, ice_heat_capacity))
end

const HCP = HeatCapacityParameters
cp_v(p::HCP)    = p.water_vapor_heat_capacity
cp_l(p::HCP)    = p.liquid_water_heat_capacity
cp_i(p::HCP)    = p.ice_heat_capacity
kappa_d(p::HCP) = p.dry_air_adiabatic_exponent

struct PhaseTransitionParameters{FT} <: AbstractThermodynamicsParameters{FT}
    reference_vaporization_enthalpy :: FT
    reference_sublimation_enthalpy  :: FT
    reference_temperature           :: FT
    triple_point_temperature        :: FT
    triple_point_pressure           :: FT
    water_freezing_temperature      :: FT
    ice_nucleation_temperature      :: FT
end

function PhaseTransitionParameters(FT = Float64;
                                   reference_vaporization_enthalpy = 2500800,
                                   reference_sublimation_enthalpy = 2834400,
                                   reference_temperature = 273.16,
                                   triple_point_temperature = 273.16,
                                   triple_point_pressure = 611.657,
                                   water_freezing_temperature = 273.16,
                                   ice_nucleation_temperature = 233)

   return PhaseTransitionParameters(convert(FT, reference_vaporization_enthalpy),
                                    convert(FT, reference_sublimation_enthalpy),
                                    convert(FT, reference_temperature),
                                    convert(FT, triple_point_temperature),
                                    convert(FT, triple_point_pressure),
                                    convert(FT, water_freezing_temperature),
                                    convert(FT, ice_nucleation_temperature))
end

const PTP = PhaseTransitionParameters
LH_v0(p::PTP)        = p.reference_vaporization_enthalpy
LH_s0(p::PTP)        = p.reference_sublimation_enthalpy
T_freeze(p::PTP)     = p.water_freezing_temperature
T_triple(p::PTP)     = p.triple_point_temperature
T_icenuc(p::PTP)     = p.ice_nucleation_temperature
press_triple(p::PTP) = p.triple_point_pressure
T_0(p::PTP)          = p.reference_temperature

struct HierarchicalThermodynamicsParameters{FT} <: AbstractThermodynamicsParameters{FT}
    constitutive       :: ConstitutiveParameters{FT}
    phase_transitions  :: PhaseTransitionParameters{FT}
    heat_capacity      :: HeatCapacityParameters{FT}
end

function HierarchicalThermodynamicsParameters(FT = Float64;
                                              constitutive = ConstitutiveParameters(FT),
                                              phase_transitions = PhaseTransitionParameters(FT),
                                              heat_capacity = HeatCapacityParameters(FT))

    return HierarchicalThermodynamicsParameters(constitutive, phase_transitions, heat_capacity)
end

const HTP = HierarchicalThermodynamicsParameters

gas_constant(p::HTP)   = gas_constant(p.constitutive)
molmass_dryair(p::HTP) = molmass_dryair(p.constitutive)
molmass_water(p::HTP)  = molmass_water(p.constitutive)
kappa_d(p::HTP)        = kappa_d(p.heat_capacity)
LH_v0(p::HTP)          = LH_v0(p.phase_transitions)
LH_s0(p::HTP)          = LH_s0(p.phase_transitions)
cp_v(p::HTP)           = cp_v(p.heat_capacity)
cp_l(p::HTP)           = cp_l(p.heat_capacity)
cp_i(p::HTP)           = cp_i(p.heat_capacity)
T_freeze(p::HTP)       = T_freeze(p.phase_transitions)
T_triple(p::HTP)       = T_triple(p.phase_transitions)
T_icenuc(p::HTP)       = T_icenuc(p.phase_transitions)
press_triple(p::HTP)   = press_triple(p.phase_transitions)
T_0(p::HTP)            = T_0(p.phase_transitions)

default_universal_function_parameters(FT=Float64) = BusingerParams{FT}(Pr_0 = convert(FT, 0.74),
                                                                       a_m  = convert(FT, 4.7),
                                                                       a_h  = convert(FT, 4.7),
                                                                       ζ_a  = convert(FT, 2.5),
                                                                       γ    = convert(FT, 4.42))

