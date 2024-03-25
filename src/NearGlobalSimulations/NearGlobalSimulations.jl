module NearGlobalSimulations

using Oceananigans
using Oceananigans.Units

using Oceananigans.Operators: Δzᵃᵃᶜ, ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Architectures: arch_array
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis
# using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
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

using ClimaOcean: u_bottom_drag, v_bottom_drag, u_immersed_bottom_drag, v_immersed_bottom_drag

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

using Oceananigans.Operators
using Oceananigans.TurbulenceClosures
using Oceananigans.Operators: Δx, Δy, ℑxyz

@inline Δ²(i, j, k, grid, lx, ly, lz)  = (1 / (1 / Δx(i, j, k, grid, lx, ly, lz)^2 + 1 / Δy(i, j, k, grid, lx, ly, lz)^2))
@inline grid_dependent_biharmonic_viscosity(i, j, k, grid, lx, ly, lz, clock, fields, τ) =
    Δ²(i, j, k, grid, lx, ly, lz)^2 / τ

@inline function leith_dynamic_biharmonic_viscosity(i, j, k, grid, lx, ly, lz, clock, fields, p)
    location = (lx, ly, lz)
    from_∂xζ = (Center(), Face(), Center())
    from_∂yζ = (Face(), Center(), Center())
    from_∂xδ = (Face(), Center(), Center())
    from_∂yδ = (Center(), Face(), Center())

    ∂xζ = ℑxyz(i, j, k, grid, from_∂xζ, location, ∂xᶜᶠᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂yζ = ℑxyz(i, j, k, grid, from_∂yζ, location, ∂yᶠᶜᶜ, ζ₃ᶠᶠᶜ, fields.u, fields.v)
    ∂xδ = ℑxyz(i, j, k, grid, from_∂xδ, location, ∂xᶠᶜᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)
    ∂yδ = ℑxyz(i, j, k, grid, from_∂yδ, location, ∂yᶜᶠᶜ, div_xyᶜᶜᶜ, fields.u, fields.v)

    dynamic_visc = sqrt( p.Cζ * (∂xζ^2 + ∂yζ^2) + p.Cδ * (∂xδ^2 + ∂yδ^2) )

    return Δ²(i, j, k, grid, lx, ly, lz)^2.5 * dynamic_visc
end

geometric_viscosity(formulation, timescale) = ScalarBiharmonicDiffusivity(formulation, ν = grid_dependent_biharmonic_viscosity, 
                                                                          discrete_form = true, 
                                                                          parameters = timescale)

function leith_viscosity(formulation; C_vort = 3.0, C_div = 3.0)

    Cζ = (C_vort / π)^6 / 8
    Cδ = (C_div  / π)^6 / 8

    return ScalarBiharmonicDiffusivity(formulation, ν = leith_dynamic_biharmonic_viscosity,
                                       discrete_form = true,
                                       parameters = (; Cζ, Cδ))
end

NoBiogeochemistry(args...; kwargs...) = nothing

include("one_degree_global_simulation.jl")
include("quarter_degree_global_simulation.jl")

end # module
