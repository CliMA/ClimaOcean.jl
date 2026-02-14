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

using Oceananigans.Fields: interpolate!
using Oceananigans.Grids: cpu_face_constructor_x,
                          cpu_face_constructor_y,
                          cpu_face_constructor_z,
                          topology

# Should we move this to grids??
construct_grid(::Type{<:RectilinearGrid}, arch, size, extent, topology) =
    RectilinearGrid(arch; size, x = extent[1], y = extent[2], z = extent[2], topology)

construct_grid(::Type{<:LatitudeLongitudeGrid}, arch, size, extent, topology) =
    LatitudeLongitudeGrid(arch; size, longitude = extent[1], latitude = extent[2], z = extent[3], topology)

include("diffuse_tracers.jl")

end # module
