module NearGlobalSimulations

using Oceananigans
using Oceananigans.Units

using Oceananigans.Operators: Δzᵃᵃᶜ, ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
using Oceananigans.Coriolis: WetCellEnstrophyConservingScheme
using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity, FluxTapering

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity, MixingLength

using DataDeps
using Statistics
using JLD2
using Printf
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using CUDA: @allowscalar

import Oceananigans.Utils: prettysummary
import ..VerticalGrids

struct PiecewiseConstantVerticalDiffusivity <: Function
    z_transition :: Float64
    top_diffusivity :: Float64
    bottom_diffusivity :: Float64
end

@inline function (pcd::PiecewiseConstantVerticalDiffusivity)(x, y, z, t)
    zᵗ = pcd.z_transition
    κᵗ = pcd.top_diffusivity
    κᵇ = pcd.bottom_diffusivity
    return ifelse(z > zᵗ, κᵗ, κᵇ)
end

prettysummary(pcd::PiecewiseConstantVerticalDiffusivity) =
    string("PiecewiseConstantVerticalDiffusivity(",
           pcd.z_transition, ", ",
           pcd.top_diffusivity, ", ",
           pcd.bottom_diffusivity, ")")

Base.summary(pcd::PiecewiseConstantVerticalDiffusivity) = prettysummary(pcd)

const thirty_days = 30days

@inline current_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
@inline next_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
@inline cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / thirty_days, 1) * (u₂ - u₁)

@inline function surface_temperature_relaxation(i, j, grid, clock, fields, p)
    time = clock.time

    n₁ = current_time_index(time, p.Nmonths)
    n₂ = next_time_index(time, p.Nmonths)

    @inbounds begin
        T★₁ = p.T★[i, j, n₁]
        T★₂ = p.T★[i, j, n₂]
        Q★₁ = p.Q★[i, j, n₁]
        Q★₂ = p.Q★[i, j, n₂]
        T_surface = fields.T[i, j, grid.Nz]
    end

    T★ = cyclic_interpolate(T★₁, T★₂, time)
    Q★ = cyclic_interpolate(Q★₁, Q★₂, time)

    return Q★ + p.λ * (T_surface - T★)
end

@inline function surface_salinity_relaxation(i, j, grid, clock, fields, p)
    time = clock.time

    n₁ = current_time_index(time, p.Nmonths)
    n₂ = next_time_index(time, p.Nmonths)

    @inbounds begin
        S★₁ = p.S★[i, j, n₁]
        S★₂ = p.S★[i, j, n₂]
        F★₁ = p.F★[i, j, n₁]
        F★₂ = p.F★[i, j, n₂]
        S_surface = fields.S[i, j, grid.Nz]
    end

    S★ = cyclic_interpolate(S★₁, S★₂, time)
    F★ = cyclic_interpolate(F★₁, F★₂, time)

    return - F★ + p.λ * (S_surface - S★)
end

@inline function surface_wind_stress(i, j, grid, clock, fields, p)
    time = clock.time
    n₁ = current_time_index(time, p.Nmonths)
    n₂ = next_time_index(time, p.Nmonths)

    @inbounds begin
        τ₁ = p.τ[i, j, n₁]
        τ₂ = p.τ[i, j, n₂]
    end

    return cyclic_interpolate(τ₁, τ₂, time)
end

@inline ϕ²(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

include("one_degree_global_simulation.jl")

end # module
