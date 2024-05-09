import Base: -
import Statistics

#####
##### Struct that represents a 3-tuple of momentum, heat, and water vapor
#####

struct SimilarityScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

-(a::SimilarityScales, b::SimilarityScales) = SimilarityScales(a.momentum - b.momentum, 
                                                                      a.temperature - b.temperature,
                                                                      a.water_vapor - b.water_vapor)

Statistics.norm(a::SimilarityScales) = norm(a.momentum) + norm(a.temperature) + norm(a.water_vapor)

# Implementation of stability functions that follow Edson et al (2013)
# We can swap them out easily if we define new types `NewStability` with a method `(f::NewStability)(ζ)`
struct MomentumStabilityFunction end
struct ScalarStabilityFunction end
struct InitialMomentumStabilityFunction end

function default_stability_functions(FT = Float64)

    # Edson et al. (2013)
    ψu = MomentumStabilityFunction()
    ψc = ScalarStabilityFunction()

    return SimilarityScales(ψu, ψc, ψc)
end

@inline function (ψ::InitialMomentumStabilityFunction)(ζ)
    FT = eltype(ζ)

    # Parameters
    p₁ = convert(FT, 0.35)

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, p₁ * ζ⁺)

    ψ_stable = - ζ⁺ - 3 / 4 * (ζ⁺ - 5 / p₁) * exp(-dζ) - 3 / 4 * 5 / p₁
    
    fₘ = sqrt(sqrt(1 - 18 * ζ⁻))
    ψ_unstable_1 = log((1 + fₘ)^2 * (1 + fₘ^2) / 8) - 2 * atan(fₘ) + π / 2;

    fₘ = cbrt(1 - 10 * ζ⁻)
    ψ_unstable_2 = 3 / 2 * log((1 + fₘ + fₘ^2) / 3) - sqrt(3) * atan((1 + 2fₘ) / sqrt(3)) + π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

@inline function (ψ::MomentumStabilityFunction)(ζ)
    FT = eltype(ζ)

    # Parameters
    p₁ = convert(FT, 0.35)
    p₂ = convert(FT, 0.7)
    p₃ = convert(FT, 10.15)

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, p₁ * ζ⁺)

    ψ_stable = - p₂ * ζ⁺ - 3 / 4 * (ζ⁺ - 5 / p₁) * exp(-dζ) - 3 / 4 * 5 / p₁
    
    fₘ = sqrt(sqrt(1 - 15 * ζ⁻))
    ψ_unstable_1 = log((1 + fₘ)^2 * (1 + fₘ^2) / 8) - 2 * atan(fₘ) + π / 2;

    fₘ = cbrt(1 - p₃ * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₘ + fₘ^2) / 3) - sqrt(3) * atan((1 + 2fₘ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

@inline function (ψ::ScalarStabilityFunction)(ζ)
    FT = eltype(ζ)

    # Parameters
    p₁ = convert(FT, 0.35)
    p₂ = convert(FT, 14.28)
    p₃ = convert(FT, 8.525)
    p₄ = convert(FT, 34.15)

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, p₁ * ζ⁺)

    ψ_stable = - (4 * ζ⁺ / 3)^(3 / 2) - 2 / 3 * (ζ⁺ - p₂) * exp(-dζ) - p₃
    
    fₕ = sqrt(1 - 15 * ζ⁻)
    ψ_unstable_1 = 2 * log((1 + fₕ) / 2) 

    fₕ = cbrt(1 - p₄ * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₕ + fₕ^2) / 3) - sqrt(3) * atan((1 + 2fₕ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end
