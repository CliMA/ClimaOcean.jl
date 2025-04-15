module InitialConditions

export initialize!

using Oceananigans
using Oceananigans: location
using Oceananigans.BoundaryConditions
using Oceananigans.Fields: OneField
using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior, CenterField
using Oceananigans.Architectures: architecture, device, GPU, child_architecture

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll
using JLD2

# Implementation of 3-dimensional regridding
# TODO: move all the following to Oceananigans! 

using Oceananigans.Fields: regrid!, interpolate!
using Oceananigans.Grids: cpu_face_constructor_x, 
                          cpu_face_constructor_y, 
                          cpu_face_constructor_z,
                          topology

include("diffuse_tracers.jl")

end # module

