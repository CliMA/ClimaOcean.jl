module Bathymetry

export regrid_bathymetry
export OceanBasinMask
export AtlanticOceanMask, IndianOceanMask, SouthernOceanMask, PacificOceanMask
export label_ocean_basins
export Barrier

using Downloads
using ImageMorphology
using KernelAbstractions: @kernel, @index
using Oceananigans
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: DistributedGrid, reconstruct_global_grid, all_reduce
using Oceananigans.Fields: interpolate!
using Oceananigans.Grids: x_domain, y_domain, topology, AbstractGrid
using Oceananigans.Utils: launch!
using OffsetArrays
using NCDatasets
using Printf
using Scratch

using ..DataWrangling: Metadatum, native_grid, metadata_path, download_dataset
using ..DataWrangling.ETOPO: ETOPO2022

include("label_ocean_basins.jl")
include("regrid_bathymetry.jl")
include("ocean_basin_mask.jl")

end # module
