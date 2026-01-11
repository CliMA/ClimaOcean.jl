using Oceananigans.Grids: λnode, φnode, on_architecture
using ImageMorphology

#####
##### OceanBasinMask struct
#####

"""
    OceanBasinMask{M, G}

A mask identifying cells belonging to a specific ocean basin.
The mask is stored as a 2D Field{Center, Center, Nothing} with values:
- 1.0: cell belongs to the basin
- 0.0: cell does not belong to the basin
"""
struct OceanBasinMask{M, G, S}
    mask :: M
    grid :: G
    seed_points :: S
end

Base.summary(obm::OceanBasinMask) = "OceanBasinMask"

function Base.show(io::IO, obm::OceanBasinMask)
    print(io, summary(obm), " on ", summary(obm.grid))
end

# Forward getindex to mask
Base.getindex(obm::OceanBasinMask, i, j, k) = obm.mask[i, j, k]

#####
##### Connected component utilities
#####

"""
    find_label_at_point(labels, grid, λs, φs)

Find the connected component label at the given longitude/latitude seed point.
Returns the label value, or 0 if the point is on land or outside the domain.
Checks in a circular cap of radius `radius` degrees around the seed point
"""
function find_label_at_point(labels, grid, λs, φs; radius = 4)
    Nx, Ny, _ = size(grid)

    # Find grid cell containing the seed point
    for j in 1:Ny
        for i in 1:Nx
            λ = λnode(i, j, 1, grid, Center(), Center(), Center())
            φ = φnode(i, j, 1, grid, Center(), Center(), Center())

            # Check if this cell contains the seed point (within a certain radius)
            # TODO: generalize this heuristic?
            if (λ - λs)^2 + (φ - φs)^2 < radius # 2 degree circular cap (enough for a basin?)
                return labels[i, j]
            end
        end
    end

    return 0  # Seed point not found
end

"""
    create_basin_mask_from_label(grid, labels, basin_label;
                                 south_boundary = nothing,
                                 north_boundary = nothing,
                                 west_boundary  = nothing,
                                 east_boundary  = nothing)

Create a mask field for all cells with the given label, optionally restricted to a latitude range 
or a longitude range.
"""
function create_basin_mask_from_label(grid, labels, basin_label;
                                      south_boundary = nothing,
                                      north_boundary = nothing,
                                      west_boundary  = nothing,
                                      east_boundary  = nothing)
    arch = architecture(grid)
    Nx, Ny, _ = size(grid)
    FT = eltype(grid)

    mask = Field{Center, Center, Nothing}(grid)

    launch!(CPU(), grid, :xy, _compute_basin_mask!, 
            mask, grid, labels, basin_label, 
            south_boundary, north_boundary,
            west_boundary, east_boundary)

    fill_halo_regions!(mask)

    return mask
end

@kernel function _compute_basin_mask!(mask, grid, labels, basin_label, south_boundary, north_boundary, west_boundary, east_boundary)
    i, j = @index(Global, NTuple)

    # Check longitude and latitude bounds if specified
    λ = λnode(i, j, 1, grid, Center(), Center(), Center())
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())

    outside_south_bounds = !isnothing(south_boundary) && (φ < south_boundary)
    outside_north_bounds = !isnothing(north_boundary) && (φ > north_boundary)

    outside_west_bounds = !isnothing(west_boundary) && (λ < west_boundary) 
    outside_east_bounds = !isnothing(east_boundary) && (λ > east_boundary)

    outside_bounds = outside_south_bounds | outside_north_bounds | outside_west_bounds | outside_east_bounds

    correct_basin = @inbounds labels[i, j] == basin_label

    @inbounds mask[i, j, 1] = !outside_bounds & correct_basin
end

#####
##### Some usefull Basin seeds and barriers
#####

# TODO: add constructors for pacific ocean and arctic ocean

const ATLANTIC_OCEAN_BARRIERS = [
    Barrier(-180.0, 180.0, -56.0, -54.0),   # Disconnect from Southern Ocean
    Barrier(20.0, -60.0, -30.0),  # Cape Agulhas (meridional barrier)
]

const INDIAN_OCEAN_BARRIERS = [
    Barrier(-180.0, 180.0, -56.0, -54.0),   # Disconnect from Southern Ocean
    Barrier(141.0, -60.0, -3.0),    # Indonesian side (meridional)
    Barrier(20.0,  -60.0, -30.0),  # Cape Agulhas (meridional barrier)
    Barrier(105.0, 141.0, -4.0, -3.0),  # Cape Agulhas (meridional barrier)
]

const SOUTHERN_OCEAN_BARRIERS = [
    LatitudeBand(-35.0, -33.0),
]

# Seed points for Atlantic Ocean (definitely in the Atlantic)
const ATLANTIC_SEED_POINTS = [
    (-30.0, 0.0),    # Central equatorial Atlantic
    (-40.0, 30.0),   # North Atlantic
    (-25.0, -20.0),  # South Atlantic
    # Same values but from 0 to 360
    (-30.0 + 360, 0.0),    # Central equatorial Atlantic
    (-40.0 + 360, 30.0),   # North Atlantic
    (-25.0 + 360, -20.0),  # South Atlantic
]

# Seed points for Indian Ocean
const INDIAN_SEED_POINTS = [
    (70.0, -10.0),   # Central Indian Ocean
    (60.0, 10.0),    # Arabian Sea region
    (90.0, -20.0),   # Eastern Indian Ocean
]

# Seed points for Southern Ocean
const SOUTHERN_SEED_POINTS = [
    (0.0, -60.0),     # South Atlantic sector
    (90.0, -60.0),    # Indian Ocean sector
    (180.0, -60.0),   # Pacific sector (date line)
    (-90.0, -60.0),   # South Pacific sector
]

##### 
##### OceanBasinMask
#####

"""
    OceanBasinMask(grid;
                   south_boundary = -34.0,
                   north_boundary = 65.0,
                   west_boundary = nothing,
                   east_boundary = nothing,
                   seed_points = ATLANTIC_SEED_POINTS,
                   barriers = nothing)

Create a mask identifying cells in the Ocean basin using connected component analysis.
The algorithm uses labels restricted by south, north, east, and west boundaries with a fill
algorithm that retrieves the basin corresponding to labels associated with known `seed_points`.
This ensures the mask exactly follows coastlines.

Arguments
=========
- `grid`: The ocean grid (typically an ImmersedBoundaryGrid)

Keyword Arguments
=================
- `south_boundary`: Southern latitude limit. Default: -34.0, at Cape of Good Hope
- `north_boundary`: Northern latitude limit. Default: 65.0
- `west_boundary`: Western longitude limit. Default: nothing
- `east_boundary`: Eastern longitude limit. Default: nothing
- `seed_points`: Known (λ, φ) points of the ocean basin to retrieve. Default: ATLANTIC_SEED_POINTS
- `barriers`: Collection of barriers to apply before labeling. Barriers temporarily
              separate connected ocean basins (e.g., separating Atlantic from Pacific).
              Use predefined barriers like `ATLANTIC_OCEAN_BARRIERS` for common basins.

Returns
=======
An `OceanBasinMask` with a 2D mask field.
"""
function OceanBasinMask(grid;
                        south_boundary = nothing,
                        north_boundary = nothing,
                        west_boundary = nothing,
                        east_boundary = nothing,
                        seed_points = [(0, 0)],
                        barriers = nothing)

    # Compute connected component labels for all ocean cells
    # Barriers are applied to temporarily separate connected basins
    labels = label_ocean_basins(grid; barriers)

    # Find the Basin label using seed points
    basin_label = 0
    for (λs, φs) in seed_points
        label = find_label_at_point(labels, grid, λs, φs)
        if label > 0
            basin_label = label
            break
        end
    end

    if basin_label == 0
        @warn "Could not find the Ocean basin in grid. Returning empty mask."
        mask = Field{Center, Center, Nothing}(grid)
        return OceanBasinMask(mask, grid, seed_points)
    end

    # Create mask from label with latitude bounds
    mask = create_basin_mask_from_label(grid, labels, basin_label;
                                        south_boundary,
                                        north_boundary,
                                        east_boundary,
                                        west_boundary)

    return OceanBasinMask(mask, grid, seed_points)
end

#####
##### Convenience functions for specific ocean basins
#####

"""
    AtlanticOceanMask(grid; kw...)

Create a mask for the Atlantic Ocean with predefined barriers and seed points.
Default boundaries: south=-34.0 (Cape of Good Hope), north=65.0
"""
function AtlanticOceanMask(grid;
                           south_boundary = -50.0,
                           north_boundary = 75.0,
                           barriers = ATLANTIC_OCEAN_BARRIERS,
                           seed_points = ATLANTIC_SEED_POINTS,
                           kw...)
    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end

"""
    IndianOceanMask(grid; kw...)

Create a mask for the Indian Ocean with predefined barriers and seed points.
Default boundaries: south=-60.0 (Southern Ocean boundary), north=30.0
"""
function IndianOceanMask(grid;
                         south_boundary = -60.0,
                         north_boundary = 30.0,
                         barriers = INDIAN_OCEAN_BARRIERS,
                         seed_points = INDIAN_SEED_POINTS,
                         kw...)
    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end

"""
    SouthernOceanMask(grid; kw...)

Create a mask for the Southern Ocean with predefined barriers and seed points.
Default boundaries: south=-90.0, north=-35.0
"""
function SouthernOceanMask(grid;
                           south_boundary = -90.0,
                           north_boundary = -35.0,
                           barriers = SOUTHERN_OCEAN_BARRIERS,
                           seed_points = SOUTHERN_SEED_POINTS,
                           kw...)
    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end
