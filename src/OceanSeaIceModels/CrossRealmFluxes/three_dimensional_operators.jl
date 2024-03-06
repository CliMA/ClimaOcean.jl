using Oceananigans.Operators
using Oceananigans.Grids: peripheral_node

@inline ı(i, j, k, grid, f::Function, args...) = f(i, j, k, grid, args...)
@inline ı(i, j, k, grid, ϕ)                    = ϕ[i, j, k]

# Defining Interpolation operators for the immersed boundaries
@inline conditional_ℑx_f(LY, LZ, i, j, k, grid, t) = ifelse(peripheral_node(i, j, k, grid, Center(), LY, LZ), ı(i-1, j, k, grid, t), ifelse(peripheral_node(i-1, j, k, grid, Center(), LY, LZ), ı(i, j, k, grid, t), ℑxᶠᵃᵃ(i, j, k, grid, t)))
@inline conditional_ℑy_f(LX, LZ, i, j, k, grid, t) = ifelse(peripheral_node(i, j, k, grid, LX, Center(), LZ), ı(i, j-1, k, grid, t), ifelse(peripheral_node(i, j-1, k, grid, LX, Center(), LZ), ı(i, j, k, grid, t), ℑyᵃᶠᵃ(i, j, k, grid, t)))
@inline conditional_ℑx_c(LY, LZ, i, j, k, grid, t) = ifelse(peripheral_node(i, j, k, grid, Face(), LY, LZ), ı(i+1, j, k, grid, t), ifelse(peripheral_node(i+1, j, k, grid, Face(), LY, LZ), ı(i, j, k, grid, t), ℑxᶜᵃᵃ(i, j, k, grid, t)))
@inline conditional_ℑy_c(LX, LZ, i, j, k, grid, t) = ifelse(peripheral_node(i, j, k, grid, LX, Face(), LZ), ı(i, j+1, k, grid, t), ifelse(peripheral_node(i, j+1, k, grid, LX, Face(), LZ), ı(i, j, k, grid, t), ℑyᵃᶜᵃ(i, j, k, grid, t)))

ℑxᶠᶜᶜ(i, j, k, grid, t) = conditional_ℑx_f(Center(), Center(), i, j, k, grid, t)
ℑyᶜᶠᶜ(i, j, k, grid, t) = conditional_ℑy_f(Center(), Center(), i, j, k, grid, t)

ℑxᶜᶜᶜ(i, j, k, grid, t) = conditional_ℑx_c(Center(), Center(), i, j, k, grid, t)
ℑyᶜᶜᶜ(i, j, k, grid, t) = conditional_ℑy_c(Center(), Center(), i, j, k, grid, t)