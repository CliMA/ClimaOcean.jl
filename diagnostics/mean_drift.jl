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

include("mixed_layer.jl")

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

Kmean  = zeros(Float32, 60, length(T.times))
Tmean  = zeros(Float32, 60, length(T.times))
Smean  = zeros(Float32, 60, length(T.times))

lu = location(u[1])
lv = location(v[1])

MetricField(loc, grid, metric) = compute!(Field(Oceananigans.AbstractOperations.GridMetricOperation(loc, metric, grid)))
VolumeField(grid, loc=(Center, Center, Center)) = MetricField(loc, grid, Oceananigans.AbstractOperations.volume)

volC  = VolumeField(grid, (Center, Center, Center));
volCF = VolumeField(grid, (Face, Center, Center));
volCF = VolumeField(grid, (Center, Face, Center));

condition_field!(volC)
condition_field!(volCF)
condition_field!(volFC)

for t in 1:length(T.times)
    @info "averaging time $t of $(length(T.times))"
    @info "day number $(T.times[t] / 86400)"
    Tmean[:, t] .= conditioned_mean(T[t], volC, dims = (1, 2))[1, 1, :]
    Smean[:, t] .= conditioned_mean(S[t], volC, dims = (1, 2))[1, 1, :]

    u2 = conditioned_mean(u[t].^2, volFC, grid, dims = (1, 2))[1, 1, :]
    v2 = conditioned_mean(v[t].^2, volCF, grid, dims = (1, 2))[1, 1, :]

    Kmean[:, t] .= 0.5 .* (u2 .+ v2)

    GC.gc()
end

jldsave("drifts.jld2"; Tmean, Smean, Kmean)

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
