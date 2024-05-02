import Base
import Statistics

#####
##### Struct that represents a 3-tuple of momentum, heat, and water vapor
#####

struct SimilarityScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

Base.(-)(a::SimilarityScales, b::SimilarityScales) = SimilarityScales(a.momentum - b.momentum, 
                                                                      a.temperature - b.temperature,
                                                                      a.water_vapor - b.water_vapor)

Statistic.norm(a::SimilarityScales) = norm(a.momentum) + norm(a.temperature) + norm(a.water_vapor)

struct MomentumStabilityFunction end
struct ScalarStabilityFunction end
struct InitialMomentumStabilityFunction end

function default_stability_functions(FT = Float64)

    # Computed from Edisen et al. (2013)
    ψu = MomentumStabilityFunction()
    ψc = ScalarStabilityFunction()

    return SimilarityScales(ψu, ψc, ψc)
end

@inline function (ψ::InitialMomentumStabilityFunction)(ζ)
    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, 0.35 * ζ⁺)

    ψ_stable = - ζ⁺ - 3 / 4 * (ζ⁺ - 5 / 0.35) * exp(-dζ) - 3 / 4 * 5 / 0.35
    
    fₘ = sqrt(sqrt(1 - 18 * ζ⁻))
    ψ_unstable_1 = log((1 + fₘ)^2 * (1 + fₘ^2) / 8) - 2 * atan(fₘ) + π / 2;

    fₘ = cbrt(1 - 10 * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₘ + fₘ^2) / 3) - sqrt(3) * atan((1 + 2fₘ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

@inline function (ψ::MomentumStabilityFunction)(ζ)

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, 0.35 * ζ⁺)

    ψ_stable = - 0.7 * ζ⁺ - 3 / 4 * (ζ⁺ - 5 / 0.35) * exp(-dζ) - 3 / 4 * 5 / 0.35
    
    fₘ = sqrt(sqrt(1 - 15 * ζ⁻))
    ψ_unstable_1 = log((1 + fₘ)^2 * (1 + fₘ^2) / 8) - 2 * atan(fₘ) + π / 2;

    fₘ = cbrt(1 - 10.15 * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₘ + fₘ^2) / 3) - sqrt(3) * atan((1 + 2fₘ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end

@inline function (ψ::ScalarStabilityFunction)(ζ)

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, 0.35 * ζ⁺)

    ψ_stable = - (4 * ζ⁺ / 3)^(3 / 2) - 2 / 3 * (ζ⁺ - 14.28) * exp(-dζ) - 8.525
    
    fₕ = sqrt(1 - 15 * ζ⁻)
    ψ_unstable_1 = 2 * log((1 + fₕ) / 2) 

    fₕ = cbrt(1 - 34.15 * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₕ + fₕ^2) / 3) - sqrt(3) * atan((1 + 2fₕ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end
