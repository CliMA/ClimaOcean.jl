module InitialConditions

export initialize!

using Oceananigans
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

# Should we move this to grids??
construct_grid(::Type{<:RectilinearGrid}, arch, size, extent, topology) = 
    RectilinearGrid(arch; size, x = extent[1], y = extent[2], z = extent[2], topology)

construct_grid(::Type{<:LatitudeLongitudeGrid}, arch, size, extent, topology) = 
    LatitudeLongitudeGrid(arch; size, longitude = extent[1], latitude = extent[2], z = extent[3], topology)

function three_dimensional_regrid!(a, b)
    target_grid = a.grid isa ImmersedBoundaryGrid ? a.grid.underlying_grid : a.grid
    source_grid = b.grid isa ImmersedBoundaryGrid ? b.grid.underlying_grid : b.grid 

    topo = topology(target_grid)
    arch = architecture(target_grid)
    arch = child_architecture(arch)
    
    target_y = yt = cpu_face_constructor_y(target_grid)
    target_z = zt = cpu_face_constructor_z(target_grid)

    target_size = Nt = size(target_grid)

    source_x = xs = cpu_face_constructor_x(source_grid)
    source_y = ys = cpu_face_constructor_y(source_grid)

    source_size = Ns = size(source_grid)

    # Start by regridding in z
    @debug "Regridding in z"
    zgrid   = construct_grid(typeof(target_grid), arch, (Ns[1], Ns[2], Nt[3]), (xs, ys, zt), topo)
    field_z = Field(location(b), zgrid)
    regrid!(field_z, zgrid, source_grid, b)

    # regrid in y 
    @debug "Regridding in y"
    ygrid   = construct_grid(typeof(target_grid), arch, (Ns[1], Nt[2], Nt[3]), (xs, yt, zt), topo)
    field_y = Field(location(b), ygrid);
    regrid!(field_y, ygrid, zgrid, field_z);

    # Finally regrid in x
    @debug "Regridding in x"
    regrid!(a, target_grid, ygrid, field_y)

    return a
end

import Oceananigans.Fields: interpolate!
using Oceananigans.Fields: _interpolate!, AbstractField
using Oceananigans.Architectures: child_architecture, architecture
using Oceananigans.Utils: launch!
using Oceananigans.BoundaryConditions
    
"""
    interpolate!(to_field::Field, from_field::AbstractField)

Interpolate `from_field` `to_field` and then fill the halo regions of `to_field`.
"""
function interpolate!(to_field::Field, from_field::AbstractField)
    to_grid   = to_field.grid
    from_grid = from_field.grid

    to_arch   = child_architecture(architecture(to_field))
    from_arch = child_architecture(architecture(from_field))
    if !isnothing(from_arch) && to_arch != from_arch
        msg = "Cannot interpolate! because from_field is on $from_arch while to_field is on $to_arch."
        throw(ArgumentError(msg))
    end

    # Make locations
    from_location = Tuple(L() for L in location(from_field))
    to_location   = Tuple(L() for L in location(to_field))

    launch!(to_arch, to_grid, size(to_field),
            _interpolate!, to_field, to_grid, to_location,
            from_field, from_grid, from_location)

    fill_halo_regions!(to_field)

    return nothing
end

function three_dimensional_interpolate!(a, b)
    target_grid = a.grid isa ImmersedBoundaryGrid ? a.grid.underlying_grid : a.grid
    source_grid = b.grid isa ImmersedBoundaryGrid ? b.grid.underlying_grid : b.grid 

    topo = topology(source_grid)
    arch = architecture(target_grid)
    arch = child_architecture(arch)
    
    target_y = yt = cpu_face_constructor_y(target_grid)
    target_z = zt = cpu_face_constructor_z(target_grid)

    target_size = Nt = size(target_grid)

    source_x = xs = cpu_face_constructor_x(source_grid)
    source_y = ys = cpu_face_constructor_y(source_grid)

    source_size = Ns = size(source_grid)

    # Start by regridding in z
    @debug "Interpolating in z"
    zgrid   = construct_grid(typeof(target_grid), arch, (Ns[1], Ns[2], Nt[3]), (xs, ys, zt), topo)
    field_z = Field(location(b), zgrid)
    interpolate!(field_z, b)

    # regrid in y 
    @debug "Interpolating in y"
    ygrid   = construct_grid(typeof(target_grid), arch, (Ns[1], Nt[2], Nt[3]), (xs, yt, zt), topo)
    field_y = Field(location(b), ygrid);
    interpolate!(field_y, field_z)

    # Finally regrid in x
    @debug "Interpolating in x"
    interpolate!(a, field_y)

    return a
end

include("diffuse_tracers.jl")

end # module

