module InterpolateData

using PyCall
using JLD2
using NetCDF
using Oceananigans.Fields: interpolate, set!
using Oceananigans.Architectures: device, arch_array, device_event
using Oceananigans.Units
using Oceananigans.Fields: interpolate
using Statistics: dot
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.Grids: architecture

export interpolate_bathymetry_from_file, remove_connected_regions, write_bathymetry_to_file!

const ABOVE_SEA_LEVEL = 100

include("utils.jl")
include("interpolate_bathymetry.jl")
# include("interpolate_initial_condition.jl")
# include("interpolate_fluxes.jl")

const scikitimage = PyNULL()

function __init__()
    copy!(scikitimage, pyimport_conda("skimage.measure", "scikit-image"))
end

end