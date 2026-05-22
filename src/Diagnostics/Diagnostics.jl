module Diagnostics

export MixedLayerDepthField, MixedLayerDepthOperand
export MeridionalStreamfunction, compute_streamfunction
export compute_amoc, compute_broken_isolatitudes, BrokenIsoLatitude

using Oceananigans
using Oceananigans.Architectures: architecture
using Oceananigans.Models: buoyancy_operation
using Oceananigans.Grids: new_data, inactive_cell, znode
using Oceananigans.BoundaryConditions: FieldBoundaryConditions, fill_halo_regions!
using Oceananigans.Fields: FieldStatus
using Oceananigans.Utils: launch!
using KernelAbstractions: @index, @kernel
using Oceananigans.Operators: ζ₃ᶠᶠᶜ, ℑxᶜᵃᵃ, ℑyᵃᶜᵃ
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Grids: on_architecture
using Oceananigans.Architectures: child_architecture
using Statistics: mean

import Oceananigans.Fields: compute!

include("mixed_layer_depth.jl")

end # module
