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

function default_stability_functions(FT = Float64)

    # Edson et al. (2013)
    ψu = MomentumStabilityFunction(FT)
    ψc = ScalarStabilityFunction(FT)

    return SimilarityScales(ψu, ψc, ψc)
end

"""
    MomentumStabilityFunction{FT}

A struct representing the momentum stability function detailed in Edson et al (2013).
The formulation hinges on the definition of three different functions:
one for stable atmospheric conditions ``(ζ > 0)``, named ``ψₛ`` and two for unstable conditions,
named ``ψᵤ₁`` and ``ψᵤ₂``.
These stability functions are obtained by regression to experimental data.

The stability parameter for stable atmospheric conditions is defined as
```math
dζ = min(ζmax, Aˢζ)
ψₛ = - (Bˢ ζ + Cˢ ( ζ - Dˢ ) ) exp( - dζ) - Cˢ Dˢ 
```

While the stability parameter for unstable atmospheric conditions is calculated
as a function of the two individual stability functions as follows
```math
fᵤ₁ = √√(1 - Aᵘζ)
ψᵤ₁ = Bᵘ / 2 ⋅ log((1 + fᵤ₁ + fᵤ₁² + fᵤ₁³) / Bᵘ) - √Bᵘ atan(fᵤ₁) - Cᵘ

fᵤ₂ = ∛(1 - Dᵘζ)
ψᵤ₂ = Eᵘ / 2 ⋅ log((1 + fᵤ₂ + fᵤ₂²) / Eᵘ) - √Eᵘ atan( (1 + 2fᵤ₂) / √Eᵘ) + Fᵘ

f  = ζ² / (1 + ζ²)
ψᵤ = (1 - f) ψᵤ₁ + f ψᵤ₂  
```

The superscripts ``ˢ`` and ``ᵘ`` indicate if the parameter applies to the 
stability function for _stable_ or _unstable_ atmospheric conditions, respectively.
"""
@kwdef struct MomentumStabilityFunction{FT}
    ζmax :: FT = 50.0
    Aˢ   :: FT = 0.35
    Bˢ   :: FT = 0.7
    Cˢ   :: FT = 0.75
    Dˢ   :: FT = 5/0.35
    Aᵘ   :: FT = 15.0
    Bᵘ   :: FT = 2.0
    Cᵘ   :: FT = π/2
    Dᵘ   :: FT = 10.15
    Eᵘ   :: FT = 3.0
    Fᵘ   :: FT = π / sqrt(3)
end

"""
    ScalarStabilityFunction{FT}

A struct representing the scalar stability function detailed in Edson et al (2013).
The formulation hinges on the definition of three different functions:
one for stable atmospheric conditions ``(ζ > 0)``, named ``ψₛ`` and two for unstable conditions,
named ``ψᵤ₁`` and ``ψᵤ₂``.
These stability functions are obtained by regression to experimental data.

The stability parameter for stable atmospheric conditions is defined as
```math
dζ = min(ζmax, Aˢζ)
ψₛ = - (1 + Bˢ ζ) ^ Cₛ - Bˢ ( ζ - Dˢ ) * exp( - dζ) - Eˢ
```

While the stability parameter for unstable atmospheric conditions is calculated
as a function of the two individual stability functions as follows
```math
fᵤ₁ = √(1 - Aᵘζ)
ψᵤ₁ = Bᵘ ⋅ log((1 + fᵤ₁) / Bᵘ) + Cᵤ

fᵤ₂ = ∛(1 - Dᵘζ)
ψᵤ₂ = Eᵘ / 2 ⋅ log((1 + fᵤ₂ + fᵤ₂²) / Eᵘ) - √Eᵘ atan( (1 + 2fᵤ₂) / √Eᵘ) + Fᵘ

f  = ζ² / (1 + ζ²)
ψᵤ = (1 - f) ψᵤ₁ + f ψᵤ₂  
```

The superscripts ``ˢ`` and ``ᵘ`` indicate if the parameter applies to the 
stability function for _stable_ or _unstable_ atmospheric conditions, respectively.
"""
@kwdef struct ScalarStabilityFunction{FT}
    ζmax :: FT = 50.0
    Aˢ   :: FT = 0.35
    Bˢ   :: FT = 0.7
    Cˢ   :: FT = 0.75
    Dˢ   :: FT = 5/0.35
    Aᵘ   :: FT = 15.0
    Bᵘ   :: FT = 4.0
    Cᵘ   :: FT = 0.0
    Dᵘ   :: FT = 10.15
    Eᵘ   :: FT = 3.0
    Fᵘ   :: FT = π / sqrt(3)
end

@inline function (ψ::MomentumStabilityFunction)(ζ)
    ζmax = ψ.ζmax
    Aˢ   = ψ.Aˢ  
    Bˢ   = ψ.Bˢ  
    Cˢ   = ψ.Cˢ  
    Dˢ   = ψ.Dˢ  
    Aᵘ   = ψ.Aᵘ  
    Bᵘ   = ψ.Bᵘ  
    Cᵘ   = ψ.Cᵘ  
    Dᵘ   = ψ.Dᵘ  
    Eᵘ   = ψ.Eᵘ  
    Fᵘ   = ψ.Fᵘ  

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(ζmax, Aˢ * ζ⁺)

    # Stability parameter for _stable_ atmospheric conditions
    ψₛ = - (Bˢ * ζ⁺ + Cˢ * (ζ⁺ - Dˢ)) * exp(- dζ) - Cˢ * Dˢ
        
    # Stability parameter for _unstable_ atmospheric conditions
    fᵤ₁ = sqrt(sqrt(1 - Aᵘ * ζ⁻))
    ψᵤ₁ = Bᵘ * log((1 + fᵤ₁) / Bᵤ) + log((1 + fᵤ₁^2) / Bᵘ) - Bᵘ * atan(fᵤ₁) + Cᵘ
        
    fᵤ₂ = cqrt(1 - Dᵘ * ζ⁻)
    ψᵤ₂ = Eᵘ / 2 * log((1 + fᵤ₂ + fᵤ₂^2) / Eᵘ) - sqrt(Eᵘ) * atan( (1 + 2fᵤ₂) / sqrt(Eᵘ)) + Fᵘ
        
    f  = ζ⁻^2 / (1 + ζ⁻^2)
    ψᵤ = (1 - f) * ψᵤ₁ + f * ψᵤ₂  
        
    return ifelse(ζ < 0, ψᵤ, ψₛ)
end

@inline function (ψ::ScalarStabilityFunction)(ζ)
    # Parameters
    p₁ = ψ.p₁
    p₂ = ψ.p₂
    p₃ = ψ.p₃
    p₄ = ψ.p₄
    p₅ = ψ.p₅

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(50, p₁ * ζ⁺)

    # stability function for stable atmospheric conditions 
    ψ_stable = - (4 * ζ⁺ / 3)^(3 / 2) - 2 / 3 * (ζ⁺ - p₂) * exp(-dζ) - p₃
    
    fₕ = sqrt(1 - p₄ * ζ⁻)
    ψ_unstable_1 = 2 * log((1 + fₕ) / 2) 

    fₕ = cbrt(1 - p₅ * ζ⁻)
    ψ_unstable_2 = 1.5 * log((1 + fₕ + fₕ^2) / 3) - sqrt(3) * atan((1 + 2fₕ) / sqrt(3))+ π / sqrt(3)
    
    f⁻ = ζ⁻^2 / (1 + ζ⁻^2)
    
    # stability function for unstable atmospheric conditions
    ψ_unstable = (1 - f⁻) * ψ_unstable_1 + f⁻ * ψ_unstable_2

    return ifelse(ζ < 0, ψ_unstable, ψ_stable)
end
