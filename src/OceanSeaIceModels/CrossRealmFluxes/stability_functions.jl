import Statistics
import Base: -

#####
##### Struct that represents a 3-tuple of momentum, heat, and water vapor
#####

struct SimilarityScales{U, T, Q}
    momentum :: U
    temperature :: T
    water_vapor :: Q
end

function -(a::SimilarityScales, b::SimilarityScales)
    Δu = a.momentum - b.momentum
    Δθ = a.temperature - b.temperature
    Δq = a.water_vapor - b.water_vapor
    return SimilarityScales(Δu, Δθ, Δq)
end

Statistics.norm(a::SimilarityScales) = norm(a.momentum) + norm(a.temperature) + norm(a.water_vapor)

# Edson et al. (2013)
function edson_stability_functions(FT = Float64)
    ψu = MomentumStabilityFunction()
    ψc = ScalarStabilityFunction()
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
    ψᵤ₁ = Bᵘ * log((1 + fᵤ₁) / Bᵘ) + log((1 + fᵤ₁^2) / Bᵘ) - Bᵘ * atan(fᵤ₁) + Cᵘ
        
    fᵤ₂ = cbrt(1 - Dᵘ * ζ⁻)
    ψᵤ₂ = Eᵘ / 2 * log((1 + fᵤ₂ + fᵤ₂^2) / Eᵘ) - sqrt(Eᵘ) * atan( (1 + 2fᵤ₂) / sqrt(Eᵘ)) + Fᵘ
        
    f  = ζ⁻^2 / (1 + ζ⁻^2)
    ψᵤ = (1 - f) * ψᵤ₁ + f * ψᵤ₂  
        
    return ifelse(ζ < 0, ψᵤ, ψₛ)
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
    Bˢ   :: FT = 2/3
    Cˢ   :: FT = 3/2
    Dˢ   :: FT = 14.28
    Eˢ   :: FT = 8.525
    Aᵘ   :: FT = 15.0
    Bᵘ   :: FT = 2.0
    Cᵘ   :: FT = 0.0
    Dᵘ   :: FT = 34.15
    Eᵘ   :: FT = 3.0
    Fᵘ   :: FT = π / sqrt(3)
end

@inline function (ψ::ScalarStabilityFunction)(ζ)
    ζmax = ψ.ζmax
    Aˢ   = ψ.Aˢ  
    Bˢ   = ψ.Bˢ  
    Cˢ   = ψ.Cˢ  
    Dˢ   = ψ.Dˢ  
    Eˢ   = ψ.Eˢ  
    Aᵘ   = ψ.Aᵘ  
    Bᵘ   = ψ.Bᵘ  
    Cᵘ   = ψ.Cᵘ  
    Dᵘ   = ψ.Dᵘ  
    Eᵘ   = ψ.Eᵘ  
    Fᵘ   = ψ.Fᵘ  

    ζ⁻ = min(zero(ζ), ζ)
    ζ⁺ = max(zero(ζ), ζ)
    dζ = min(ζmax, Aˢ * ζ⁺)

    # stability function for stable atmospheric conditions 
    ψₛ = - (1 + Bˢ * ζ⁺)^Cˢ - Bˢ * (ζ⁺ - Dˢ) * exp(-dζ) - Eˢ

    # Stability parameter for _unstable_ atmospheric conditions
    fᵤ₁ = sqrt(1 - Aᵘ * ζ⁻)
    ψᵤ₁ = Bᵘ * log((1 + fᵤ₁) / Bᵘ) + Cᵘ

    fᵤ₂ = cbrt(1 - Dᵘ * ζ⁻)
    ψᵤ₂ = Eᵘ / 2 * log((1 + fᵤ₂ + fᵤ₂^2) / Eᵘ) - sqrt(Eᵘ) * atan( (1 + 2fᵤ₂) / sqrt(Eᵘ)) + Fᵘ

    f  = ζ⁻^2 / (1 + ζ⁻^2)
    ψᵤ = (1 - f) * ψᵤ₁ + f * ψᵤ₂  

    return ifelse(ζ < 0, ψᵤ, ψₛ)
end
