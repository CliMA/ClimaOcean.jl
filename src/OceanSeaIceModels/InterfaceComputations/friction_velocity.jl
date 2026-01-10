"""
    MomentumBasedFrictionVelocity

A friction velocity formulation that computes the friction velocity from momentum stresses.

The friction velocity is computed as:
```math
u_* = \\sqrt{\\frac{|\\boldsymbol{\\tau}|}{\\rho_o}}
```
where τ is the magnitude of the momentum stress vector and ρₒ is the ocean reference density.

Example
=======

```jldoctest
using ClimaOcean.OceanSeaIceModels: MomentumBasedFrictionVelocity

fv = MomentumBasedFrictionVelocity()

# output
MomentumBasedFrictionVelocity (computed from momentum stresses)
```
"""
struct MomentumBasedFrictionVelocity end

@inline ϕ²(i, j, k, grid, ϕ)        = @inbounds ϕ[i, j, k]^2
@inline τᶜᶜᶜ(i, j, k, grid, τx, τy) = @inbounds sqrt(ℑxᶜᵃᵃ(i, j, k, grid, ϕ², τx) + ℑyᵃᶜᵃ(i, j, k, grid, ϕ², τy))

Base.summary(::MomentumBasedFrictionVelocity) = "MomentumBasedFrictionVelocity"

function Base.show(io::IO, ::MomentumBasedFrictionVelocity)
    print(io, "MomentumBasedFrictionVelocity (computed from momentum stresses)")
end

"""
    get_friction_velocity(u★, i, j, args...)

Return the friction velocity at grid point `(i, j)`.

For a constant friction velocity (`u★::Number`), returns the value directly.
For `MomentumBasedFrictionVelocity`, computes ``u_* = \\sqrt{|\\tau| / \\rho_o}`` from momentum stresses.
"""
@inline get_friction_velocity(u★::Number, i, j, args...) = u★
@inline get_friction_velocity(::MomentumBasedFrictionVelocity, i, j, τx, τy, ρₒ) = sqrt(τᶜᶜᶜ(i, j, 1, τx, τy) / ρ₀)
