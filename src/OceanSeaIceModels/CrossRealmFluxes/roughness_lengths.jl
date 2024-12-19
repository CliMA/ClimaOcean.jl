struct MomentumRoughnessLength{FT, V}
    gravitational_acceleration :: FT
    air_kinematic_viscosity :: V
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
    maximum_roughness_length :: FT
end

struct ScalarRoughnessLength{FT, V, R}
    air_kinematic_viscosity :: V
    reynolds_number_scaling_function :: R
    maximum_roughness_length :: FT
end

"""
    ScalarRoughnessLength([FT=Float64];
                          air_kinematic_viscosity = temperature_dependent_viscosity,
                          reynolds_number_scaling_function = empirical_scaling_function,
                          maximum_roughness_length = 1.6e-4)

Constructs a `ScalarRoughnessLength` object that represents the scalar roughness length
that regulates the exchange of heat and water vapor between the ocean and the atmosphere.   

Keyword Arguments
==================

- `air_kinematic_viscosity::Function`: The function to compute the air kinematic viscosity.
- `reynolds_number_scaling_function::Function`: The function to compute the Reynolds number scaling factor.
- `maximum_roughness_length::Float`: The maximum roughness length value. Defaults to `1.6e-4`.
"""
function ScalarRoughnessLength(FT=Float64;
                               air_kinematic_viscosity = TemperatureDependentAirViscosity(FT),
                               reynolds_number_scaling_function = ReynoldsScalingFunction(FT),
                               maximum_roughness_length = 1.6e-4) # Values from COARE3.6

    return ScalarRoughnessLength(air_kinematic_viscosity,
                                 reynolds_number_scaling_function,
                                 convert(FT, maximum_roughness_length))
end

"""
    MomentumRoughnessLength([FT=Float64];
                            gravitational_acceleration = default_gravitational_acceleration,
                            maximum_roughness_length = 1.0,
                            air_kinematic_viscosity = TemperatureDependentAirViscosity(FT),
                            gravity_wave_parameter = 0.011,
                            laminar_parameter = 0.11)

Construct a `MomentumRoughnessLength` object that represents the momentum roughness length that
regulates the exchange of momentum, heat, and water vapor between the ocean and the atmosphere.

Keyword Arguments
=================

- `gravitational_acceleration`: The gravitational acceleration. Default: `default_gravitational_acceleration`.
- `maximum_roughness_length`: The maximum roughness length. Default: 1.0.
- `air_kinematic_viscosity`: The air kinematic viscosity. Default: `TemperatureDependentAirViscosity(FT)`.
- `gravity_wave_parameter`: The wave parameter. Default: 0.011.
- `laminar_parameter`: The laminar parameter. Default: 0.11.
"""
function MomentumRoughnessLength(FT=Float64;
                                 gravitational_acceleration = default_gravitational_acceleration,
                                 maximum_roughness_length = 1.0, # An estimate?
                                 air_kinematic_viscosity = TemperatureDependentAirViscosity(FT),
                                 gravity_wave_parameter = 0.011,
                                 laminar_parameter = 0.11)

    return MomentumRoughnessLength(convert(FT, gravitational_acceleration),
                                          air_kinematic_viscosity,
                                          convert(FT, gravity_wave_parameter),
                                          convert(FT, laminar_parameter),
                                          convert(FT, maximum_roughness_length))
end

function default_roughness_lengths(FT=Float64)
    momentum    = MomentumRoughnessLength(FT)
    temperature = ScalarRoughnessLength(FT)
    water_vapor = ScalarRoughnessLength(FT)
    return SimilarityScales(momentum, temperature, water_vapor)
end

# Temperature-dependent viscosity law
struct TemperatureDependentAirViscosity{FT}
    C₀ :: FT
    C₁ :: FT
    C₂ :: FT
    C₃ :: FT
end

"""
    TemperatureDependentAirViscosity([FT = Float64;
                                      C₀ = 1.326e-5,
                                      C₁ = C₀ * 6.542e-3,
                                      C₂ = C₀ * 8.301e-6,
                                      C₃ = - C₀ * 4.84e-9])

Constructs a `TemperatureDependentAirViscosity` object that calculates the kinematic
viscosity of air as 
```math
C₀ + C₁ T + C₂ T^2 + C₃ T^3.
```
"""
function TemperatureDependentAirViscosity(FT = Float64;
                                          C₀ = 1.326e-5,
                                          C₁ = C₀ * 6.542e-3,
                                          C₂ = C₀ * 8.301e-6,
                                          C₃ = - C₀ * 4.84e-9)

    return TemperatureDependentAirViscosity(convert(FT, C₀),
                                            convert(FT, C₁),
                                            convert(FT, C₂),
                                            convert(FT, C₃))
end   

""" Calculate the air viscosity based on the temperature θ in Celsius. """
@inline function (ν::TemperatureDependentAirViscosity)(θ)
    FT = eltype(ν.C₀)
    T  = convert(FT, θ - celsius_to_kelvin)
    return ν.C₀ + ν.C₁ * T + ν.C₂ * T^2 + ν.C₃ * T^3
end

# Fallbacks for constant roughness length!
@inline roughness_length(ℓ, u★, args...)     = ℓ(u★, args...)
@inline roughness_length(ℓ::Number, args...) = ℓ

# Momentum roughness length should be different from scalar roughness length.
# Temperature and water vapor can be considered the same (Edson et al 2013)
@inline function roughness_length(ℓ::MomentumRoughnessLength{FT}, u★, 𝒬, ℂ) where FT
    g  = ℓ.gravitational_acceleration
    α  = ℓ.gravity_wave_parameter
    β  = ℓ.laminar_parameter
    ℓm = ℓ.maximum_roughness_length

    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬)
    ν  = ℓ.air_kinematic_viscosity(θ₀)

    # We need to prevent `Inf` that pops up when `u★ == 0`.
    # For this reason, if `u★ == 0` we prescribe the roughness length to be
    # equal to a `maximum` roughness length
    ℓᴿ = ifelse(u★ == 0, ℓm, β * ν / u★) 
    
    return min(α * u★^2 / g + ℓᴿ, ℓm)
end

struct ReynoldsScalingFunction{FT}
    A :: FT
    b :: FT
end

"""
    ReynoldsScalingFunction(FT = Float64; A = 5.85e-5, b = 0.72)

Empirical fit of the scalar roughness length with roughness Reynolds number `R★ = u★ ℓu / ν`.
Edson et al. (2013), equation (28).
```math
    ℓs = A / R★ ^ b
```
"""
ReynoldsScalingFunction(FT = Float64; A = 5.85e-5, b = 0.72) = 
    ReynoldsScalingFunction(convert(FT, A), convert(FT, b))

@inline (s::ReynoldsScalingFunction)(R★, args...) = ifelse(R★ == 0, convert(eltype(R★), 0), s.A / R★ ^ s.b)

# Edson 2013 formulation of scalar roughness length
@inline function roughness_length(ℓ::ScalarRoughnessLength{FT}, ℓu, u★, 𝒬, ℂ) where FT
    ℓm = ℓ.maximum_roughness_length
    
    scaling_function = ℓ.reynolds_number_scaling_function

    θ₀ = AtmosphericThermodynamics.air_temperature(ℂ, 𝒬)
    ν  = ℓ.air_kinematic_viscosity(θ₀)

    # Roughness Reynolds number
    R★ = ℓu * u★ / ν

    # implementation of scalar roughness length
    ℓq = scaling_function(R★, ℓu, u★, ν)

    return min(ℓq, ℓm) 
end
