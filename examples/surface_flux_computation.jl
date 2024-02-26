 using SurfaceFluxes
 using Thermodynamics
 using StaticArrays
 using ClimaOcean

 import CLIMAParameters

 using Thermodynamics: q_vap_saturation_from_density, partial_pressure_vapor
 using SurfaceFluxes.Parameters: SurfaceFluxesParameters
 using ClimaOcean.OceanSeaIceModels:
    default_universal_function_parameters,
    default_surface_flux_parameters

 const CP = CLIMAParameters

 function extrapolate_surface_density(params, atmos_state, surface_temperature)
     Tₛ = surface_temperature
     Tₐ = air_temperature(params, atmos_state)
     Rmₐ = gas_constant_air(params, atmos_state)
     ρₐ = air_density(params, atmos_state)
     κ = cv_m(params, atmos_state) / Rmₐ
     return ρₐ * (Tₛ / Tₐ)^κ
end

function surface_saturation_specific_humidity(params, surface_temperature, atmos_state)
    Tₛ = surface_temperature
    ρₛ = atmos_state.ρ #extrapolate_surface_density(thermo_params, atmos_state, Tₛ)

    @show p★ = saturation_vapor_pressure(params, Tₛ, Liquid())

    q = PhasePartition(thermo_params, atmos_thermo_state)
    @show q

    @show pᵥ = partial_pressure_vapor(params, atmos_state.p, q)
    @show p★ₐ = saturation_vapor_pressure(params, atmos_state.T, Liquid())
    @show q★ₐ = q_vap_saturation_from_density(params, atmos_state.T, atmos_state.ρ, p★ₐ)
    @show pᵥ / p★ₐ

    q★ = q_vap_saturation_from_density(params, Tₛ, ρₛ, p★)
    @show q★
    @show q

    return q★
end
 
FT = Float64
thermo_params = Thermodynamics.Parameters.HierarchicalThermodynamicsParameters(FT)
businger_params = default_universal_function_parameters(FT)
surface_flux_parameters = default_surface_flux_parameters(thermo_params)

#=
include(joinpath(pkgdir(SurfaceFluxes), "parameters", "create_parameters.jl"))
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
uf_type = SurfaceFluxes.UniversalFunctions.BusingerType()
param_set = create_parameters(toml_dict, uf_type)
thermo_params = SurfaceFluxes.Parameters.thermodynamics_params(param_set)
=#

h = 2.0 # height at which measurements are made, in m
surface_velocity = SVector(0.0, 0.0)
atmos_velocity = SVector(4.0, 0.0)

atmos_pressure = 101350.0
atmos_temperature = 298.15
atmos_specific_humidity = 0.03

atmos_thermo_state = Thermodynamics.PhaseEquil_pTq(thermo_params,
                                                   atmos_pressure,
                                                   atmos_temperature,
                                                   atmos_specific_humidity)


surface_temperature = Tₛ = 297.15

q₀ = 0.98
c₁ = 640380
c₂ = 5107.4
qLY = q₀ * c₁ * exp(-c₂ / Tₛ)

q★ = surface_saturation_specific_humidity(thermo_params, surface_temperature, atmos_thermo_state)
surface_thermo_state = Thermodynamics.PhaseEquil_pTq(thermo_params,
                                                     atmos_pressure,
                                                     surface_temperature,
                                                     q★)


# State at z=0, eg the "surface"
surface_dynamic_state = SurfaceFluxes.StateValues(0.0, surface_velocity, surface_thermo_state)

# State at z=h, eg the "atmosphere"
atmos_dynamic_state = SurfaceFluxes.StateValues(h, atmos_velocity, atmos_thermo_state)

momentum_roughness_length = 0.01
buoyancy_roughness_length = 0.001

values = SurfaceFluxes.ValuesOnly(atmos_dynamic_state,
                                  surface_dynamic_state,
                                  momentum_roughness_length,
                                  buoyancy_roughness_length)

conditions = SurfaceFluxes.surface_conditions(surface_flux_parameters, values)

