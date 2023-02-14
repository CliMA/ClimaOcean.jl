__precompile__()
module Initialization

using PyCall
using JLD2
using Oceananigans.Fields: interpolate, set!
using Oceananigans.Architectures: device, arch_array, device_event
using Oceananigans.Units
using Oceananigans.Fields: interpolate
using Statistics: dot
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.Grids: architecture

# an small number larger than 0, 
# to avoid finite precision errors
const ABOVE_SEA_LEVEL = 0.001

include("utils.jl")
include("interpolate_bathymetry.jl")
include("interpolate_initial_condition.jl")
include("interpolate_fluxes.jl")

const sckikitimage = PyNULL()

function __init__()
    copy!(sckikitimage, pyimport_conda("skimage.measure", "scikit-image"))
end

end # module
