
####
#### This file contains all the bug-fixes that still didn't get merged in Oceananigans.jl
####

using Oceananigans.BuoyancyModels: ∂z_b
using Oceananigans.Operators
using Oceananigans.Grids: peripheral_node, inactive_node
using Oceananigans.TurbulenceClosures: top_buoyancy_flux, getclosure, taper

import Oceananigans.TurbulenceClosures: _compute_ri_based_diffusivities!

@inline function _compute_ri_based_diffusivities!(i, j, k, diffusivities, grid, closure,
                                                  velocities, tracers, buoyancy, tracer_bcs, clock)

    # Ensure this works with "ensembles" of closures, in addition to ordinary single closures
    closure_ij = getclosure(i, j, closure)

    ν₀  = closure_ij.ν₀
    κ₀  = closure_ij.κ₀
    κᶜᵃ = closure_ij.κᶜᵃ
    Cᵃᵛ = closure_ij.Cᵃᵛ
    Ri₀ = closure_ij.Ri₀
    Riᵟ = closure_ij.Riᵟ
    tapering = closure_ij.Ri_dependent_tapering

    # Convection and entrainment
    N² = ∂z_b(i, j, k, grid, buoyancy, tracers)

    # Conditions
    convecting = N² < 0 # applies regardless of Qᵇ

    # Convective adjustment diffusivity
    κᶜᵃ = ifelse(convecting, κᶜᵃ, zero(grid))

    # Shear mixing diffusivity and viscosity
    Ri = ℑxyᶜᶜᵃ(i, j, k, grid, ℑxyᶠᶠᵃ, diffusivities.Ri)
    τ = taper(tapering, Ri, Ri₀, Riᵟ)
    
    κᶜ★ = κ₀ * τ
    κᵘ★ = ν₀ * τ

    # Previous diffusivities
    κᶜ = diffusivities.κᶜ
    κᵘ = diffusivities.κᵘ

    # New diffusivities
    κᶜ⁺ = κᶜ★ + κᶜᵃ
    κᵘ⁺ = κᵘ★

    # Set to zero on periphery and NaN within inactive region
    on_periphery = peripheral_node(i, j, k, grid, Center(), Center(), Face())
    within_inactive = inactive_node(i, j, k, grid, Center(), Center(), Face())
    κᶜ⁺ = ifelse(on_periphery, zero(grid), ifelse(within_inactive, NaN, κᶜ⁺))
    κᵘ⁺ = ifelse(on_periphery, zero(grid), ifelse(within_inactive, NaN, κᵘ⁺))

    # Update by averaging in time
    @inbounds κᶜ[i, j, k] = (Cᵃᵛ * κᶜ[i, j, k] + κᶜ⁺) / (1 + Cᵃᵛ)
    @inbounds κᵘ[i, j, k] = (Cᵃᵛ * κᵘ[i, j, k] + κᵘ⁺) / (1 + Cᵃᵛ)
    
    return nothing
end