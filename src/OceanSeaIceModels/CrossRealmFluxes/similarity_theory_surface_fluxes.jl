using Oceananigans.Utils: prettysummary

using SurfaceFluxes.Parameters: SurfaceFluxesParameters, AbstractSurfaceFluxesParameters
using SurfaceFluxes.UniversalFunctions: BusingerParams

using ..PrescribedAtmospheres: PrescribedAtmosphereThermodynamicsParameters

import Thermodynamics as AtmosphericThermodynamics

import SurfaceFluxes.Parameters:
    thermodynamics_params,
    uf_params,
    von_karman_const,
    universal_func_type,
    grav

#####
##### Bulk turbulent fluxes based on similarity theory
#####

struct SimilarityTheoryTurbulentFluxes{FT, ΔU, UF, TP} <: AbstractSurfaceFluxesParameters
    gravitational_acceleration :: FT
    von_karman_constant :: FT
    bulk_velocity_scale :: ΔU
    universal_function :: UF
    thermodynamics_parameters :: TP
end

const STTF = SimilarityTheoryTurbulentFluxes
thermodynamics_params(fluxes::STTF) = fluxes.thermodynamics_parameters
universal_func_type(fluxes::STTF)   = universal_func_type(typeof(fluxes.universal_function))
uf_params(fluxes::STTF)             = fluxes.universal_function
von_karman_const(fluxes::STTF)      = fluxes.von_karman_constant
grav(fluxes::STTF)                  = fluxes.gravitational_acceleration

Base.summary(::SimilarityTheoryTurbulentFluxes{FT}) where FT = "SimilarityTheoryTurbulentFluxes{$FT}"

function Base.show(io::IO, fluxes::SimilarityTheoryTurbulentFluxes)
    print(io, summary(fluxes), '\n',
          "├── gravitational_acceleration: ", prettysummary(fluxes.gravitational_acceleration), '\n',
          "├── von_karman_constant: ",        prettysummary(fluxes.von_karman_constant), '\n',
          "├── bulk_velocity_scale: ",        summary(fluxes.bulk_velocity_scale), '\n',
          "├── universal_function: ",         summary(fluxes.universal_function), '\n',
          "└── thermodynamics_parameters: ",  summary(fluxes.thermodynamics_parameters))
end

const PATP = PrescribedAtmosphereThermodynamicsParameters

function SimilarityTheoryTurbulentFluxes(FT = Float64;
                                         gravitational_acceleration = 9.80665,
                                         bulk_velocity_scale = nothing,
                                         von_karman_constant = 0.4,
                                         universal_function = default_universal_function_parameters(FT),
                                         thermodynamics_parameters = PATP(FT))

    return SimilarityTheoryTurbulentFluxes(gravitational_acceleration,
                                           von_karman_constant,
                                           bulk_velocity_scale,
                                           universal_function,
                                           thermodynamics_parameters)
end

# See SurfaceFluxes.jl for other parameter set options.
default_universal_function_parameters(FT=Float64) = BusingerParams{FT}(Pr_0 = convert(FT, 0.74),
                                                                       a_m  = convert(FT, 4.7),
                                                                       a_h  = convert(FT, 4.7),
                                                                       ζ_a  = convert(FT, 2.5),
                                                                       γ    = convert(FT, 4.42))

function surface_saturation_specific_humidity(params, surface_temperature, atmos_state, surface_type)
    Tₛ = surface_temperature
    ρₛ = atmos_state.ρ # should we extrapolate to obtain this?
    p★ = AtmosphericThermodynamics.saturation_vapor_pressure(params, Tₛ, surface_type)
    q★ = AtmosphericThermodynamics.q_vap_saturation_from_density(params, Tₛ, ρₛ, p★)
    q★ *= 0.98 # TODO: understand this better...
    return q★
end

