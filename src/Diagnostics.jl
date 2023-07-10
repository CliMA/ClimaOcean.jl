module Diagnostics

export MixedLayerDepthField

using Oceananigans
using Oceananigans.BuoyancyModels: buoyancy
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Architectures: device, architecture
using Oceananigans.Utils: launch!
using Oceananigans.Grids: Center, Face, inactive_node, znode
using Oceananigans.Operators: Δzᶜᶜᶠ

using KernelAbstractions: @index, @kernel
using KernelAbstractions.Extras.LoopInfo: @unroll

import Oceananigans.Fields: compute!

const c = Center()
const f = Face()

@inline z_bottom(i, j, grid) = znode(c, c, f, i, j, 1, grid)
@inline z_bottom(i, j, grid::ImmersedBoundaryGrid) = @inbounds grid.immersed_boundary.bottom_height[i, j]

@inline bottom(i, j, grid) = znode(c, c, f, i, j, 1, grid)
#####
##### MixedLayerDepthField
#####

@kernel function compute_mld!(h, grid, b, Δb)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    z_surface = znode(c, c, f, i, j, Nz+1, grid)
    b_surface = @inbounds b[i, j, Nz] # buoyancy at surface
    
    # Start downwards look
    z_ij = z_surface
    found_mixed_layer_depth = false
    
    @unroll for k in Nz-1 : -1 : 1 # scroll from point just below surface
        b⁺ = @inbounds b[i, j, k+1]
        bᵏ = @inbounds b[i, j, k]

        # If buoyancy decreases downwards, both are > 0
        Δb⁺ = b_surface - b⁺
        Δbᵏ = b_surface - bᵏ

        zᵏ = znode(c, c, c, i, j, k, grid)
        Δz⁺ = Δzᶜᶜᶠ(i, j, k+1, grid)

        # Assuming buoyancy decreases downwards
        just_below_mixed_layer = (Δb⁺ < Δb) & (Δbᵏ >= Δb)
        just_below_mixed_layer *= !found_mixed_layer_depth

        # Linearly interpolate to find mixed layer height
        new_z_ij = zᵏ + (Δb - Δbᵏ) / (Δb⁺ - Δbᵏ) * Δz⁺

        # Replace z_ij if we found a new mixed layer depth
        z_ij = ifelse(!found_mixed_layer_depth | just_below_mixed_layer, new_z_ij, z_ij)
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    h_ij = z_surface - z_ij
    @inbounds h[i, j, 1] = ifelse(inactive_node(i, j, Nz, grid, c, c, c), zero(grid), h_ij)
end

struct MixedLayerDepthOperand{B, FT}
    buoyancy_operation :: B
    mixed_layer_buoyancy_differential :: FT
end

Base.summary(op::MixedLayerDepthOperand) = "MixedLayerDepthOperand"

const MixedLayerDepthField = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

"""
    MixedLayerDepthField(grid, tracers, buoyancy_model; Δb = 3e-4, kw...)

Return a reduced `Field{Center, Center, Nothing}` that represents
mixed layer depth for `model`, based on a buoyancy differential criterion.
The mixed layer depth is defined as the depth ``h`` for which

```math
b(z=0) - b(z=-h) = Δb
```

This criterion is solved by integrating downwards and linearly interpolating to find `h`,
assuming that ``b`` decreases with depth.

Keyword arguments
=================

* `Δb`: Buoyancy differential used to calculate mixed layer depth
* `field_kw`: Keyword arguments passed to `Field`.

Example
=======

```julia
h = MixedLayerDepth(model)
compute!(h) # compute mixed layer depth
```
"""
function MixedLayerDepthField(grid, tracers, buoyancy_model; Δb = 3e-4, kw...)
    b_op = buoyancy(buoyancy_model, grid, tracers)
    operand = MixedLayerDepthOperand(b_op, Δb)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

MixedLayerDepthField(model; kw...) = MixedLayerDepthField(model.grid,
                                                          model.tracers,
                                                          model.buoyancy; kw...)

function compute!(h::MixedLayerDepthField, time=nothing)
    arch = architecture(h)
    b = h.operand.buoyancy_operation
    Δb = h.operand.mixed_layer_buoyancy_differential
    launch!(arch, h.grid, :xy, compute_mld!, h, h.grid, b, Δb)
    fill_halo_regions!(h)
    return h
end

end # module
