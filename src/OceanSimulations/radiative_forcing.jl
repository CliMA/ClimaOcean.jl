using Oceananigans.Operators: ∂zᶜᶜᶜ
using Adapt

struct TwoColorRadiation{FT, J}
    first_color_fraction :: FT
    first_absorption_coefficient :: FT
    second_absorption_coefficient :: FT
    surface_flux :: J
end

Adapt.adapt_structure(to, R::TwoColorRadiation) =
    TwoColorRadiation(adapt(to, R.first_color_fraction),
                      adapt(to, R.first_absorption_coefficient),
                      adapt(to, R.second_absorption_coefficient),
                      adapt(to, R.surface_flux))

"""
    TwoColorRadiation(grid;
                      first_color_fraction,
                      first_absorption_coefficient,
                      second_absorption_coefficient)

Return `TwoColorRadiation` that computes the radiative flux divergence associated with
a two-color radiation flux that decays according to Beer's law,

```math
I(z) = ϵ₁ * I₀ * exp(κ₁ * z) + (1 - ϵ₁) * I₀ * exp(κ₂ * z)
```

where ``I₀`` is the surface flux, ``ϵ₁`` is the "first color" fraction of the total radiation,
and ``κ₁`` and ``κ₂`` are the absorption coefficients for the two colors.
"""
function TwoColorRadiation(grid;
                           first_color_fraction,
                           first_absorption_coefficient,
                           second_absorption_coefficient)
    FT = eltype(grid)
    surface_flux = Field{Center, Center, Nothing}(grid)
    return TwoColorRadiation(convert(FT, first_color_fraction),
                             convert(FT, first_absorption_coefficient),
                             convert(FT, second_absorption_coefficient),
                             surface_flux)
end

const c = Center()
const f = Face()

@inline function beers_law_radiation(i, j, k, grid, J⁰, κ)
    z = Oceananigans.Grids.znode(i, j, k, grid, c, c, f)
    return J⁰ * exp(κ * z)
end

@inline function (R::TwoColorRadiation)(i, j, k, grid, clock, fields)
    J₀ = @inbounds R.surface_flux[i, j, 1]
    κ₁ = R.first_absorption_coefficient
    κ₂ = R.second_absorption_coefficient
    ϵ₁ = R.first_color_fraction

    # Radiation flux divergences
    dJ₁dz = ∂zᶜᶜᶜ(i, j, k, grid, beers_law_radiation, J₀, κ₁)
    dJ₂dz = ∂zᶜᶜᶜ(i, j, k, grid, beers_law_radiation, J₀, κ₂)

    # Net radiation flux divergence
    return ϵ₁ * dJ₁dz + (1 - ϵ₁) * dJ₂dz
end
