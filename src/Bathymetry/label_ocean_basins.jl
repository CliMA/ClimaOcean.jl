#####
##### Barrier type for separating ocean basins
#####

"""
    Barrier{FT}

A rectangular geographic region used to separate ocean basins during labeling.
All cells within this region are temporarily marked as land.

Fields
======
- `west`: Western longitude limit (degrees)
- `east`: Eastern longitude limit (degrees)
- `south`: Southern latitude limit (degrees)
- `north`: Northern latitude limit (degrees)

Constructors
============
- `Barrier(west, east, south, north)`: Create a barrier with explicit bounds
- `Barrier(; west, east, south, north)`: Keyword argument version
"""
struct Barrier{FT}
    west  :: FT
    east  :: FT
    south :: FT
    north :: FT
end

Barrier(; west, east, south, north) = Barrier(west, east, south, north)

"""
    Barrier(longitude, south, north; width=2.0)

Create a meridional (north-south) barrier at a given longitude.
Useful for closing straits like Cape Agulhas or Indonesian passages.
"""
Barrier(longitude, south, north; width=2.0) = Barrier(longitude - width/2, longitude + width/2, south, north)

#####
##### Barrier application functions
#####

"""
    apply_barrier!(zb_data, grid, barrier::Barrier)

Mark all cells within the barrier region as land (z = 0).
"""
apply_barrier!(zb, grid, barrier::Barrier) = 
    launch!(architecture(grid), grid, :xy, _apply_barrier!, zb, grid, barrier)

apply_barrier!(zb, grid, barriers::Nothing) = zb

function apply_barrier!(zb, grid, barriers::AbstractVector)
    for barrier in barriers
        apply_barrier!(zb, grid, barrier)
    end
    return zb
end

@kernel function _apply_barrier!(zb, grid, barrier::Barrier)
    i, j = @index(Global, NTuple)

    # If the barrier spans all longitudes (360° or more), skip longitude check.
    # This handles latitudinal barriers correctly regardless of grid longitude convention.
    full_longitude_span = (barrier.east - barrier.west) >= 360

    λ = λnode(i, j, 1, grid, Center(), Center(), Center())
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())

    in_lon = full_longitude_span | (barrier.west <= λ <= barrier.east)
    in_lat = barrier.south <= φ <= barrier.north

    @inbounds zb[i, j, 1] = ifelse(in_lon & in_lat, zero(grid), zb[i, j, 1])
end

# Since the strel algorithm in `remove_major_basins` does not recognize periodic boundaries,
# before removing connected regions, we extend the longitude direction if it is periodic.
# An extension of half the domain is enough.
function maybe_extend_longitude(zb_cpu::Field, ::Periodic)
    Nx = size(zb_cpu, 1)
    nx = Nx ÷ 2

    zb_data   = zb_cpu.data[1:Nx, :, 1]
    zb_parent = zb_data.parent

    # Add information on the LHS and to the RHS
    zb_parent = vcat(zb_parent[nx:Nx, :], zb_parent, zb_parent[1:nx, :])

    # Update offsets: the original domain starts at parent index (Nx - nx + 2)
    # For offset index 1 to map there, we need offset = -(nx + 1)
    yoffsets = zb_cpu.data.offsets[2]
    xoffsets = - nx - 1

    return OffsetArray(zb_parent, xoffsets, yoffsets)
end

maybe_extend_longitude(zb_cpu::Field, tx) = interior(zb_cpu, :, :, 1)

"""
    label_ocean_basins(zb_field, TX, core_size)

Creates a matrix with a unique label for each connected basin. Useful for inpainting the bathymetry and
computing the masks for oceanic basins.

Handles periodic boundary extension internally and returns labels for the core region only.
"""
function label_ocean_basins(zb_field, TX, size)
    zb = maybe_extend_longitude(zb_field, TX()) # Outputs a 2D AbstractArray

    water = zb .< 0

    connectivity = ImageMorphology.strel(water)
    labels = ImageMorphology.label_components(connectivity)

    Nx, Ny = size[1], size[2]
    core_labels = labels[1:Nx, 1:Ny]

    # Enforce periodicity: merge labels that should be connected across
    # the periodic boundary. This handles cases where a barrier (e.g., blocking
    # the Southern Ocean) prevents the extended domain from connecting basins
    # that are actually periodic neighbors.
    enforce_periodic_labels!(core_labels, TX())

    return core_labels
end

"""
    enforce_periodic_labels!(labels, ::Periodic)

Merge labels that should be connected due to periodic boundary conditions.
For each latitude, if the westernmost and easternmost cells are both water
(non-zero labels), they are periodic neighbors and must have the same label.
"""
function enforce_periodic_labels!(labels, ::Periodic)
    Nx, Ny = size(labels)

    for j in 1:Ny
        label_west = labels[1, j]
        label_east = labels[Nx, j]

        # Both cells are water and have different labels: merge them
        if label_west != 0 && label_east != 0 && label_west != label_east
            # Replace all occurrences of label_east with label_west
            replace!(labels, label_east => label_west)
        end
    end

    return labels
end

# No-op for non-periodic domains
enforce_periodic_labels!(labels, tx) = labels

# Utilities to label ocean basins passing only the grid
function label_ocean_basins(grid::AbstractGrid; barriers=nothing)
    @warn "The grid is not immersed, there is only one ocean basin!"
    Nx, Ny, Nz = size(grid)
    return zeros(Int, Nx, Ny)
end

"""
    label_ocean_basins(grid::ImmersedBoundaryGrid; barriers=nothing)

Label connected ocean basins in an ImmersedBoundaryGrid.

Keyword Arguments
=================
- `barriers`: Collection of barriers to apply before labeling. Barriers temporarily
              mark certain cells as land, allowing separation of connected ocean basins
              (e.g., separating Atlantic from Pacific via the Southern Ocean).
"""
function label_ocean_basins(grid::ImmersedBoundaryGrid; barriers=nothing)
    
    # The labelling algorithm works only on CPUs
    cpu_grid = on_architecture(CPU(), grid)
    
    TX = topology(cpu_grid, 1)
    zb = cpu_grid.immersed_boundary.bottom_height

    # If barriers are specified, apply them to a copy of the bathymetry
    if !isnothing(barriers)
        # Create a temporary field with the modified bathymetry
        zb_modified = Field{Center, Center, Nothing}(cpu_grid)
        parent(zb_modified) .= parent(zb)
        apply_barrier!(zb, cpu_grid, barriers)
    else
        zb_modified = zb
    end

    return label_ocean_basins(zb_modified, TX, size(cpu_grid))
end
