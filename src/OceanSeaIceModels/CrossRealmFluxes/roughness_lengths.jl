struct GravityMomentumRoughnessLength{FT, V}
    gravitational_acceleration :: FT
    air_kinematic_viscosity :: V
    gravity_wave_parameter :: FT
    laminar_parameter :: FT
    maximum_roughness_length :: FT
end

struct GravityScalarRoughnessLength{FT, V, R}
    air_kinematic_viscosity :: V
    reynolds_number_scaling_function :: R
    maximum_roughness_length :: FT
end

function default_roughness_lengths(FT=Float64)
    momentum    = GravityMomentumRoughnessLength(FT)
    temperature = GravityScalarRoughnessLength(FT)
    water_vapor = GravityScalarRoughnessLength(FT)
    return SimilarityScales(momentum, temperature, water_vapor)
end

# Empirical fit of the scalar roughness length with roughness Reynolds number `R★ = u★ / ν`
# Edson et al. (2013), equation (28)
@inline empirical_scaling_function(R★ :: FT, args...) where FT = 
        ifelse(R★ == 0, FT(0), convert(FT, 5.85e-5 / R★ ^ 0.72))

# Assumes that θ comes in in Kelvin
@inline function temperature_dependent_viscosity(θ :: FT) where FT 
    T = convert(FT, θ - celsius_to_kelvin)
    ν = convert(FT, 1.326e-5 * (1 + 6.542e-3 * T + 8.301e-6 * T^2 - 4.84e-9 * T^3))
    
    return ν
end

# Fallbacks for constant roughness length!
@inline roughness_length(ℓ, u★, args...)     = ℓ(u★, args...)
@inline roughness_length(ℓ::Number, args...) = ℓ

function GravityScalarRoughnessLength(FT=Float64;
                                      air_kinematic_viscosity = temperature_dependent_viscosity,
                                      reynolds_number_scaling_function = empirical_scaling_function,
                                      maximum_roughness_length = 1.6e-4) # Values from COARE3.6

    return GravityScalarRoughnessLength(air_kinematic_viscosity,
                                        reynolds_number_scaling_function,
                                        convert(FT, maximum_roughness_length))
end

function GravityMomentumRoughnessLength(FT=Float64;
                                        gravitational_acceleration = default_gravitational_acceleration,
                                        maximum_roughness_length = 1.0, # An estimate?
                                        air_kinematic_viscosity = temperature_dependent_viscosity,
                                        gravity_wave_parameter = 0.011,
                                        laminar_parameter = 0.11)

    return GravityMomentumRoughnessLength(convert(FT, gravitational_acceleration),
                                          air_kinematic_viscosity,
                                          convert(FT, gravity_wave_parameter),
                                          convert(FT, laminar_parameter),
                                          convert(FT, maximum_roughness_length))
end

# Momentum roughness length should be different from scalar roughness length.
# Apparently temperature and water vapor can be considered the same (Edison et al 2013)
@inline function roughness_length(ℓ::GravityMomentumRoughnessLength{FT}, u★, 𝒬, ℂ) where FT
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

# This, for example is what is implemented in COARE 3.6
@inline function roughness_length(ℓ::GravityScalarRoughnessLength{FT}, ℓu, u★, 𝒬, ℂ) where FT
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
