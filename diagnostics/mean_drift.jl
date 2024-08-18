using Oceananigans
using OrthogonalSphericalShellGrids
using Oceananigans.Fields: condition_operand, ConstantField
using Oceananigans.AbstractOperations: materialize_condition!, ComputedField
using Oceananigans.AbstractOperations: compute_at!, compute_computed_field!
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Architectures: device, architecture
using Oceananigans.Utils: launch!
using Oceananigans.Models: seawater_density
using Oceananigans.Grids: Center, Face, inactive_node, znode
using Oceananigans.Operators: Δzᶜᶜᶠ, ζ₃ᶠᶠᶜ
using JLD2

using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll

import Oceananigans.Fields: compute!

using Oceananigans.Fields: OneField, condition_operand
using Oceananigans.AbstractOperations: materialize_condition!
using Oceananigans.Utils

filename = "snapshots.jld2"

function compute!(comp::ComputedField, time=nothing)
    # First compute `dependencies`:
    compute_at!(comp.operand, time)

    # Now perform the primary computation
    @apply_regionally compute_computed_field!(comp)

    return comp
end

T = FieldTimeSeries(filename, "T"; backend = OnDisk())
S = FieldTimeSeries(filename, "S"; backend = OnDisk())
u = FieldTimeSeries(filename, "u"; backend = OnDisk())
v = FieldTimeSeries(filename, "v"; backend = OnDisk())
e = FieldTimeSeries(filename, "e"; backend = OnDisk())

grid = T.grid

# Upper 700 meters
z700  = findfirst(x -> x < -700,  grid.zᵃᵃᶜ)
z2000 = findfirst(x -> x < -2000, grid.zᵃᵃᶜ)

import Oceananigans.Fields: interior

interior(u::Union{Array, CuArray, SubArray}) = u

conditioned_mean(T::Field, condition = nothing; kwargs...) = 
    conditioned_mean(interior(T), location(T), T.grid; kwargs...)

function conditioned_mean(T, loc, grid, condition = nothing; dims = Colon())
    onefield = Field(loc, grid)
    fill!(onefield, 1)
    condition_field!(onefield, condition)

    Tnew = deepcopy(T)

    condition_field!(Tnew, condition)

    mean_active = sum(interior(onefield); dims)
    Tmean       = sum(Tnew; dims)
    
    return Tmean ./ mean_active
end

function condition_field!(f, condition = nothing)
    cf = condition_operand(f, condition, 0)
    materialize_condition!(cf);
    return nothing
end

Kmean  = zeros(Float32, 60, length(T.times))
Tmean  = zeros(Float32, 60, length(T.times))
Smean  = zeros(Float32, 60, length(T.times))

lu = location(u[1])
lv = location(v[1])

for t in 1:length(T.times)
    @info "averaging time $t of $(length(T.times))"
    @info "day number $(T.times[t] / 86400)"
    Tmean[:, t] .= conditioned_mean(T[t], dims = (1, 2))[1, 1, :]
    Smean[:, t] .= conditioned_mean(S[t], dims = (1, 2))[1, 1, :]

    u2 = conditioned_mean(u[t].^2, lu, grid, dims = (1, 2))
    v2 = conditioned_mean(v[t].^2, lu, grid, dims = (1, 2))

    Kmean[:, t] .= 0.5 .* (u2 .+ v2)

    GC.gc()
end

jldsave("drifts.jld2"; Tmean, Smean, Kmean)

const c = Center()
const f = Face()

@inline z_bottom(i, j, grid) = znode(i, j, 1, grid, c, c, f)
@inline z_bottom(i, j, grid::ImmersedBoundaryGrid) = @inbounds grid.immersed_boundary.bottom_height[i, j]

@inline bottom(i, j, grid) = znode(i, j, 1, grid, c, c, f)

#####
##### MixedLayerDepthField
#####

# b can be temperature (T) or density (ρ)
@kernel function compute_mld!(h, grid, b, Δb, density_criterion)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz

    b_surface = @inbounds (b[i, j, Nz-2] + b[i, j, Nz-1] + b[i, j, Nz]) / 3

    k_start = Nz - 3
    z_ij    = znode(i, j, k_start, grid, c, c, f)

    @unroll for k in k_start : -1 : 1 # scroll from point just below surface

        b⁺ = @inbounds b[i, j, k+1]
        bᵏ = @inbounds b[i, j, k]

        # If temperature decreases downwards, both are > 0
        # If density increases downwards, both are < 0
        Δb⁺ = b_surface - b⁺
        Δbᵏ = b_surface - bᵏ

        zᵏ = znode(i, j, k, grid, c, c, c)
        Δz⁺ = Δzᶜᶜᶠ(i, j, k+1, grid)

        # Assuming temperature decreases downwards and density increases upwards
        if density_criterion # density criterion
            inside_mixed_layer = (Δb⁺ > - Δb) & (Δbᵏ > - Δb)
            just_below_mixed_layer = (Δb⁺ > - Δb) & (Δbᵏ <= - Δb)
            new_z_ij = zᵏ + (- Δb - Δbᵏ) / (Δb⁺ - Δbᵏ) * Δz⁺
        else # temperature criterion
            inside_mixed_layer = (Δb⁺ < Δb) & (Δbᵏ < Δb)
            just_below_mixed_layer = (Δb⁺ < Δb) & (Δbᵏ >= Δb)
            new_z_ij = zᵏ + (Δb - Δbᵏ) / (Δb⁺ - Δbᵏ) * Δz⁺
        end

        # Linearly interpolate to find mixed layer height

        # Replace z_ij if we found a new mixed layer depth
        replace_z = (just_below_mixed_layer | inside_mixed_layer) & !inactive_node(i, j, k, grid, c, c, c)
        z_ij = ifelse(replace_z, new_z_ij, z_ij)
        if just_below_mixed_layer
            break
        end
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    @inbounds h[i, j, 1] = ifelse(inactive_node(i, j, k_start, grid, c, c, c), zero(grid), - z_ij)
end

struct MixedLayerDepthOperand{B, FT}
    buoyancy_operation :: B
    mixed_layer_buoyancy_differential :: FT
    density_criterion :: Bool
end

Base.summary(op::MixedLayerDepthOperand) = "MixedLayerDepthOperand"

function MixedLayerDepth(grid, tracers; Δb = 0.03, density_criterion = false, kw...)
    b_op    = density_criterion ? tracers.ρ : tracers.T
    operand = MixedLayerDepthOperand(b_op, abs(Δb), density_criterion)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

const MixedLayerDepthField = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

function compute!(h::MixedLayerDepthField, time=nothing)
    arch = architecture(h)
    b    = h.operand.buoyancy_operation
    Δb   = h.operand.mixed_layer_buoyancy_differential
    crit = h.operand.density_criterion
    launch!(arch, h.grid, :xy, compute_mld!, h, h.grid, b, Δb, crit)
    return h
end

####
#### MLD
####

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

eos    = TEOS10EquationOfState()
height = ConstantField(-200)

mld = []
for t in 1:length(T.times)
    ρ_operation = seawater_density(grid, eos, T[t], S[t], height)
    ρ = compute!(Field(ρ_operation))
    mld_t = MixedLayerDepth(grid, (; ρ), density_criterion = true)
    compute!(mld_t);

    push!(mld, interior(mld_t))
end

jldsave("mixed_layer.jld2"; mld)

####
#### AMOC?
####

# MetricField(loc, grid, metric) = compute!(Field(Oceananigans.AbstractOperations.GridMetricOperation(loc, metric, grid)))
# VolumeField(grid, loc=(Center, Center, Center)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume)

# volF = VolumeField(grid, (Center, Face, Center));
# ΔxF  = MetricField((Center, Face, Center), grid, Oceananigans.AbstractOperations.Δx);

# function AMOC(v)
#     right_of_pacific(i, j, k, grid) = (j > - i + 490) | (i > 150)
#     left_of_africa(i, j, k, grid)   = (i < 500)

#     atlantic(i, j, k, grid) = right_of_pacific(i, j, k, grid) & left_of_africa(i, j, k, grid)

#     loc = location(v);
#     onefield = Field(loc, v.grid);
#     fill!(onefield, 1);
#     condition_field!(onefield, atlantic);
#     vnew = deepcopy(v);
#     condition_field!(vnew, atlantic);

#     Vatlantic = Field{Nothing, loc[2], loc[3]}(v.grid);
#     mean_active = sum(interior(onefield), dims = 1);
#     vmean = sum(interior(vnew) .* interior(ΔxF), dims = 1);

#     set!(Vatlantic,  vmean); # ./ mean_active);
#     return Vatlantic
# end

# # ψ = ∫ᶻ V dz
# @kernel function _AMOC!(ψ, grid, V)
#     j = @index(Global, Linear)

#     # At the bottom!
#     ψ[1, j, 1] = 0
#     for k in 2:grid.Nz+1
#         ψ[1, j, k] = ψ[1, j, k-1] - Δzᶜᶠᶜ(1, j, k-1, grid) * V[1, j, k-1]
#     end
# end

# interior(V)[isnan.(interior(V))] .= 0

# ψ  = Field{Nothing, Face, Face}(grid);
# ψ1 = Field{Nothing, Face, Face}(grid1);

# λ, φ, z = nodes(v)
# λ1, φ1, z1 = nodes(v1)
# zF = grid.zᵃᵃᶠ[1:101]

# _AMOC!(KernelAbstractions.CPU(), 16, grid.Ny  + 1)(ψ, grid, V);