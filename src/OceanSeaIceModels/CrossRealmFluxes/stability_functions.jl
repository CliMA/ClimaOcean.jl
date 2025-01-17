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

#=
function -(a::SimilarityScales, b::SimilarityScales)
    Δu = a.momentum - b.momentum
    Δθ = a.temperature - b.temperature
    Δq = a.water_vapor - b.water_vapor
    return SimilarityScales(Δu, Δθ, Δq)
end
=#

# Statistics.norm(a::SimilarityScales) = norm(a.momentum) + norm(a.temperature) + norm(a.water_vapor)

# Edson et al. (2013)
function edson_stability_functions(FT = Float64)
    ψu = EdsonMomentumStabilityFunctionsFunction{FT}()
    ψc = EdsonScalarStabilityFunction{FT}()
    return SimilarityScales(ψu, ψc, ψc)
end

"""
    EdsonMomentumStabilityFunctionsFunction{FT}

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
@kwdef struct EdsonMomentumStabilityFunctionsFunction{FT}
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

@inline function (ψ::EdsonMomentumStabilityFunctionsFunction)(ζ)
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
    EdsonScalarStabilityFunction{FT}

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
@kwdef struct EdsonScalarStabilityFunction{FT}
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

@inline function (ψ::EdsonScalarStabilityFunction)(ζ)
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
    ψᵤ₂ = Eᵘ / 2 * log((1 + fᵤ₂ + fᵤ₂^2) / Eᵘ) - sqrt(Eᵘ) * atan((1 + 2fᵤ₂) / sqrt(Eᵘ)) + Fᵘ

    f  = ζ⁻^2 / (1 + ζ⁻^2)
    ψᵤ = (1 - f) * ψᵤ₁ + f * ψᵤ₂  

    return ifelse(ζ < 0, ψᵤ, ψₛ)
end

#####
##### From Grachev et al 2007, for stable boundary layers
#####

@kwdef struct ShebaMomentumStabilityFunction{FT}
    a :: FT = 6.5
    b :: FT = 1.3
end

# @inline (ψ::ShebaMomentumStabilityFunction)(ζ) = 1 + ψ.a * ζ * cbrt(1 + ζ) / (ψ.b + ζ)
@inline function (Ψ::ShebaMomentumStabilityFunction)(ζ)
    a = Ψ.a
    b = Ψ.b
    ζ⁺ = max(zero(ζ), ζ)
    z = cbrt(1 + ζ⁺)
    B = cbrt((1 - b) / b)

    rt3 = sqrt(3)
    Ψ₁ = - 3 * a * (z - 1) / b
    Ψ₂ = a * B / 2b * (2 * log((z + B) / (1 + B))
                       - log((z^2 - B * z + B^2) / (1 - B + B^2))
                       + 2 * rt3 * (atan((2z - B) / (rt3 * B)) - atan((2 - B) / (rt3 * B))))

    return Ψ₁ + Ψ₂
end

@kwdef struct ShebaScalarStabilityFunction{FT}
    a :: FT = 5.0
    b :: FT = 5.0
    c :: FT = 3.0
end

@inline function (Ψ::ShebaScalarStabilityFunction)(ζ)
    a = Ψ.a
    b = Ψ.b
    c = Ψ.c
    B = sqrt(c^2 - 4)
    ζ⁺ = max(zero(ζ), ζ)

    Ψ₁ = - b/2 * log(1 + c * ζ⁺ + ζ⁺^2)
    Ψ₂ = (b * c / 2B - a / B) * (log((2ζ⁺ + c - B) / (2ζ⁺ + c + B))
                                 + log((c - B) / (c + B)))

    return Ψ₁ + Ψ₂
end

#####
##### From Paulson 1970 for unstable boundary layers
####

@kwdef struct PaulsonMomentumStabilityFunction{FT}
    a :: FT = 16.0
    b :: FT = π/2
end

@inline function (Ψ::PaulsonMomentumStabilityFunction)(ζ)
    a = Ψ.a
    b = Ψ.b
    ζ⁻ = min(zero(ζ), ζ)
    z = sqrt(sqrt((1 - a * ζ⁻)))

    Ψ₁ = 2 * log((1 + z) / 2)
    Ψ₂ = log((1 + z^2) / 2)
    Ψ₃ = - 2 * atan(z)

    return Ψ₁ + Ψ₂ + Ψ₃ + b
end

@kwdef struct PaulsonScalarStabilityFunction{FT}
    a :: FT = 16.0
end

@inline function (Ψ::PaulsonScalarStabilityFunction)(ζ)
    a = Ψ.a
    ζ⁻ = min(zero(ζ), ζ)
    z = sqrt(sqrt((1 - a * ζ⁻)))
    return 2 * log((1 + z^2) / 2)
end

struct SplitStabilityFunction{S, U}
    stable :: S
    unstable :: U
end

@inline function (Ψ::SplitStabilityFunction)(ζ)
    Ψ_stable = Ψ.stable(ζ)
    Ψ_unstable = Ψ.unstable(ζ)
    stable = ζ > 0
    return ifelse(stable, Ψ_stable, Ψ_unstable)
end

function atmosphere_sea_ice_stability_functions(FT=Float64)
    stable_momentum = PaulsonMomentumStabilityFunction{FT}()
    unstable_momentum = ShebaMomentumStabilityFunction{FT}()
    momentum = SplitStabilityFunction(stable_momentum, unstable_momentum)

    stable_scalar = PaulsonScalarStabilityFunction{FT}()
    unstable_scalar = ShebaScalarStabilityFunction{FT}()
    scalar = SplitStabilityFunction(stable_scalar, unstable_scalar)

    return SimilarityScales(momentum, scalar, scalar)
end

