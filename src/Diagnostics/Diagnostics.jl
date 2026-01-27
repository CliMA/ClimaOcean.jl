module Diagnostics

export MixedLayerDepthField, MixedLayerDepthOperand
export BrokenLineSet, IsoLatitudeBrokenLine
export LatitudinalBandTags  # Alias for BrokenLineSet
export atlantic_latitudinal_bands
export compute_streamfunction, compute_amoc_streamfunction
export band_latitudes, streamfunction_depths

using Oceananigans
using Oceananigans.Architectures: architecture
using Oceananigans.Models: buoyancy_operation
using Oceananigans.Grids: new_data, inactive_cell, znode
using Oceananigans.BoundaryConditions: FieldBoundaryConditions, fill_halo_regions!
using Oceananigans.Fields: FieldStatus
using Oceananigans.Utils: launch!
using KernelAbstractions: @index, @kernel

using ..Bathymetry: atlantic_ocean_mask

import Oceananigans.Fields: compute!

include("mixed_layer_depth.jl")
include("latitudinal_band_tagging.jl")

end # module
