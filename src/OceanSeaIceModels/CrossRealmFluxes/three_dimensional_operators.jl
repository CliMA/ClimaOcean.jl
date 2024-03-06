using Oceananigans.Operators
using Oceananigans.ImmersedBoundaries: c, f

@inline ı(i, j, k, grid, f::Function, args...) = f(i, j, k, grid, args...)
@inline ı(i, j, k, grid, ϕ)                    = ϕ[i, j, k]

# Defining Interpolation operators for the immersed boundaries
@inline conditional_ℑx_f(LY, LZ, i, j, k, grid, t) = ifelse(inactive_node(i, j, k, grid, c, LY, LZ), ı(i-1, j, k, grid, t), ifelse(inactive_node(i-1, j, k, grid, c, LY, LZ), ı(i, j, k, grid, t), ℑxᶠᵃᵃ(i, j, k, grid, t)))
@inline conditional_ℑy_f(LX, LZ, i, j, k, grid, t) = ifelse(inactive_node(i, j, k, grid, LX, c, LZ), ı(i, j-1, k, grid, t), ifelse(inactive_node(i, j-1, k, grid, LX, c, LZ), ı(i, j, k, grid, t), ℑyᵃᶠᵃ(i, j, k, grid, t)))

ℑxᶠᶜᶜ(i, j, k, grid, t) = conditional_ℑx_f(c, c, i, j, k, grid, t)
ℑyᶜᶠᶜ(i, j, k, grid, t) = conditional_ℑy_f(c, c, i, j, k, grid, t)