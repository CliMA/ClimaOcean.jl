module InitialConditions

export initialize!

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields: OneField
using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior, CenterField
using Oceananigans.Architectures: architecture, device, GPU

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll
using JLD2

# TODO: move all the following to Oceananigans!

using Oceananigans.Fields: regrid!
using Oceananigans.Grids: cpu_face_constructor_x, 
                          cpu_face_constructor_y, 
                          cpu_face_constructor_z,
                          topology

# Should we move this to grids??
construct_grid(::Type{<:RectilinearGrid}, arch, size, extent, topology) = 
    RectilinearGrid(arch; size, x = extent[1], y = extent[2], z = extent[2], topology)

construct_grid(::Type{<:LatitudeLongitudeGrid}, arch, size, extent, topology) = 
    LatitudeLongitudeGrid(arch; size, longitude = extent[1], latitude = extent[2], z = extent[3], topology)

# Extend this to regrid automatically
function three_dimensional_regrid!(a, b)
    target_grid = a.grid isa ImmersedBoundaryGrid ? a.grid.underlying_grid : a.grid
    source_grid = b.grid isa ImmersedBoundaryGrid ? b.grid.underlying_grid : b.grid 

    topo = topology(target_grid)
    arch = architecture(target_grid)
    
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
end

include("diffuse_tracers.jl")
include("initialize_model.jl")

end # module

