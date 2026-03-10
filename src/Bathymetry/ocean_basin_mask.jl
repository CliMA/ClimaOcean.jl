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
##### Additive basin masks
#####

"""
    Base.:+(mask1::OceanBasinMask, mask2::OceanBasinMask)

Combine two ocean basin masks, returning a new mask that is the union of both.
Cells that belong to either basin will have value 1.0.
"""
function Base.:+(mask1::OceanBasinMask, mask2::OceanBasinMask)
    grid = mask1.grid
    combined_mask = Field{Center, Center, Nothing}(grid)

    # Union: cell is in combined mask if it's in either mask
    parent(combined_mask) .= max.(parent(mask1.mask), parent(mask2.mask))
    fill_halo_regions!(combined_mask)

    # Combine seed points
    combined_seeds = (mask1.seed_points..., mask2.seed_points...)

    return OceanBasinMask(combined_mask, grid, combined_seeds)
end

#####
##### Connected component utilities
#####

"""
    find_label_at_point(labels, grid, λs, φs; radius = 2)

Find the connected component label at the given longitude/latitude seed point.
Returns the label value, or 0 if the point is on land or outside the domain.
Checks in a circular cap of radius `radius` degrees around the seed point
"""
function find_label_at_point(labels, grid, λs, φs; radius = 2)
    Nx, Ny, _ = size(grid)

    # Find grid cell containing the seed point
    for j in 1:Ny
        for i in 1:Nx
            λ = λnode(i, j, 1, grid, Center(), Center(), Center())
            φ = φnode(i, j, 1, grid, Center(), Center(), Center())

            λ  = convert_to_0_360(λ)

            Δλ = if isnothing(λs)
                zero(λ)
            else
                λs = convert_to_0_360(λs)
                λ - λs
            end

            # Check if this cell contains the seed point (within a certain radius)
            # TODO: generalize this heuristic?
            if Δλ^2 + (φ - φs)^2 < radius^2 # 2 degree circular cap (enough for a basin?)
                return labels[i, j]
            end
        end
    end

    return 0  # Seed point not found
end

"""
    create_basin_mask_from_label(grid, labels, basin_label)

Create a mask field for all cells with the given label, optionally restricted to a latitude range 
or a longitude range.
"""
function create_basin_mask_from_label(grid, labels, basin_label)
    arch = architecture(grid)
    Nx, Ny, _ = size(grid)
    FT = eltype(grid)

    mask = Field{Center, Center, Nothing}(grid)

    launch!(CPU(), grid, :xy, _compute_basin_mask!, 
            mask, grid, labels, basin_label)

    fill_halo_regions!(mask)

    return mask
end

@kernel function _compute_basin_mask!(mask, grid, labels, basin_label)
    i, j = @index(Global, NTuple)
    correct_basin = @inbounds labels[i, j] == basin_label
    @inbounds mask[i, j, 1] = correct_basin
end

#####
##### Some usefull Basin seeds and barriers
#####

const SOUTHERN_OCEAN_SEPARATION_BARRIER = Barrier(-180.0, 180.0, -56.0, -54.0)

const ATLANTIC_OCEAN_BARRIERS = [
    Barrier(20.0, -90.0, -30.0),      # Cape Agulhas (meridional barrier)
    Barrier(289.0, -90.0, -30.0)      # Drake passage (meridional barrier)       
]

const INDIAN_OCEAN_BARRIERS = [
    Barrier(141.0, -90.0, -3.0),        # Indonesian side (meridional)
    Barrier(20.0,  -90.0, -30.0),       # Cape Agulhas (meridional barrier)
    Barrier(105.0, 141.0, -4.0, -3.0),  # Indonesian/Asian seas (zonal barrier at 3.5ᵒ S)
]

const SOUTHERN_OCEAN_BARRIERS = [SOUTHERN_OCEAN_SEPARATION_BARRIER]

const PACIFIC_OCEAN_BARRIERS = [
    Barrier(141.0, -90.0, -3.0),        # Indonesian side (meridional)
    Barrier(20.0,  -90.0, -30.0),       # Cape Agulhas (meridional barrier)
    Barrier(105.0, 141.0, -4.0, -3.0),  # Indonesian/Asian seas (zonal barrier at 3.5ᵒ S)
]

# Seed points for Atlantic Ocean (definitely in the Atlantic)
const ATLANTIC_SEED_POINTS = [
    (-30.0, 0.0),    # Central equatorial Atlantic
    (-40.0, 30.0),   # North Atlantic
    (-25.0, -20.0),  # South Atlantic
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

# Seed points for Pacific Ocean
const PACIFIC_SEED_POINTS = [
    (180.0, 0.0),     # Central equatorial Pacific (dateline)
    (-150.0, 20.0),   # North Pacific (Hawaii region)
    (-120.0, -20.0),  # South Pacific
    # Same values but from 0 to 360
    (180.0, 0.0),             # Central equatorial Pacific
    (-150.0 + 360, 20.0),     # North Pacific
    (-120.0 + 360, -20.0),    # South Pacific
]

#####
##### OceanBasinMask
#####

add_barrier(v::AbstractVector, b::Barrier) = [v..., b]
add_barrier(v::Barrier,        b::Barrier) = [v, b]
add_barrier(::Nothing,         b::Barrier) = b

"""
    OceanBasinMask(grid;
                   south_boundary = nothing,
                   north_boundary = nothing,
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
- `south_boundary`: Southern latitude limit. Default: nothing
- `north_boundary`: Northern latitude limit. Default: nothing
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
                        seed_points = [(0, 0)],
                        barriers = nothing)

    # The computations are 2D and require serial algorithms, so
    # we perform the computation on the CPU then move the output
    # to the GPU if the initial grid was a GPU grid
    cpu_grid = Oceananigans.on_architecture(CPU(), grid)

    # Enforce north and south boundaries
    if !isnothing(south_boundary) 
        barriers = add_barrier(barriers, Barrier(nothing, nothing, -90, south_boundary))
    end

    if !isnothing(north_boundary) 
        barriers = add_barrier(barriers, Barrier(nothing, nothing, north_boundary, 90))
    end
    
    # Compute connected component labels for all ocean cells
    # Barriers are applied to temporarily separate connected basins
    labels = label_ocean_basins(cpu_grid; barriers)

    # Find the Basin label using seed points
    basin_label = 0
    for (λs, φs) in seed_points
        label = find_label_at_point(labels, cpu_grid, λs, φs)
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
    mask = create_basin_mask_from_label(cpu_grid, labels, basin_label)
    mask = Oceananigans.on_architecture(architecture(grid), mask)

    return OceanBasinMask(mask, grid, seed_points)
end

#####
##### Convenience functions for specific ocean basins
#####

"""
    atlantic_ocean_mask(grid; include_southern_ocean=false, kw...)

Create a mask for the Atlantic Ocean with predefined barriers and seed points.

Keyword Arguments
=================
- `include_southern_ocean`: If `true`, extends the Atlantic basin into the Southern Ocean
                            sector below the standard separation latitude (~55°S). Default: `false`.
- `south_boundary`: Southern latitude limit. Default: -50.0 (or -90.0 if `include_southern_ocean=true`)
- `north_boundary`: Northern latitude limit. Default: 75.0
- Other keyword arguments are passed to `OceanBasinMask`.
"""
function atlantic_ocean_mask(grid;
                             include_southern_ocean = true,
                             south_boundary = include_southern_ocean ? -90.0 : -50.0,
                             north_boundary = 65.0,
                             barriers = ATLANTIC_OCEAN_BARRIERS,
                             seed_points = ATLANTIC_SEED_POINTS,
                             kw...)

    if !include_southern_ocean
        barriers = [barriers..., SOUTHERN_OCEAN_SEPARATION_BARRIER]
    end

    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end

"""
    indian_ocean_mask(grid; include_southern_ocean=false, kw...)

Create a mask for the Indian Ocean with predefined barriers and seed points.

Keyword Arguments
=================
- `include_southern_ocean`: If `true`, extends the Indian basin into the Southern Ocean
                            sector below the standard separation latitude (~55°S). Default: `false`.
- `south_boundary`: Southern latitude limit. Default: -50.0 (or -90.0 if `include_southern_ocean=true`)
- `north_boundary`: Northern latitude limit. Default: 30.0
- Other keyword arguments are passed to `OceanBasinMask`.
"""
function indian_ocean_mask(grid;
                           include_southern_ocean = true,
                           south_boundary = include_southern_ocean ? -90.0 : -50.0,
                           north_boundary = 30.0,
                           barriers = INDIAN_OCEAN_BARRIERS,
                           seed_points = INDIAN_SEED_POINTS,
                           kw...)

    if !include_southern_ocean
        barriers = [barriers..., SOUTHERN_OCEAN_SEPARATION_BARRIER]
    end

    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end

"""
    southern_ocean_mask(grid; kw...)

Create a mask for the Southern Ocean with predefined barriers and seed points.
Default boundaries: south=-90.0, north=-35.0
"""
function southern_ocean_mask(grid;
                           south_boundary = -90.0,
                           north_boundary = -35.0,
                           barriers = SOUTHERN_OCEAN_BARRIERS,
                           seed_points = SOUTHERN_SEED_POINTS,
                           kw...)
    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end

"""
    pacific_ocean_mask(grid; include_southern_ocean=false, kw...)

Create a mask for the Pacific Ocean with predefined barriers and seed points.

Keyword Arguments
=================
- `include_southern_ocean`: If `true`, extends the Pacific basin into the Southern Ocean
                            sector below the standard separation latitude (~55°S). Default: `false`.
- `south_boundary`: Southern latitude limit. Default: -50.0 (or -90.0 if `include_southern_ocean=true`)
- `north_boundary`: Northern latitude limit. Default: 65.0
- Other keyword arguments are passed to `OceanBasinMask`.
"""
function pacific_ocean_mask(grid;
                            include_southern_ocean = true,
                            south_boundary = include_southern_ocean ? -90.0 : -50.0,
                            north_boundary = 65.0,
                            barriers = PACIFIC_OCEAN_BARRIERS,
                            seed_points = PACIFIC_SEED_POINTS,
                            kw...)

    if !include_southern_ocean
        barriers = [barriers..., SOUTHERN_OCEAN_SEPARATION_BARRIER]
    end

    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end

"""
    arctic_ocean_mask(grid; include_southern_ocean=false, kw...)

Create a mask for the Arctic Ocean with predefined barriers and seed points.

Keyword Arguments
=================
- `include_southern_ocean`: If `true`, extends the Pacific basin into the Southern Ocean
                            sector below the standard separation latitude (~55°S). Default: `false`.
- `south_boundary`: Southern latitude limit. Default: -50.0 (or -90.0 if `include_southern_ocean=true`)
- `north_boundary`: Northern latitude limit. Default: 65.0
- Other keyword arguments are passed to `OceanBasinMask`.
"""
function arctic_ocean_mask(grid;
                           include_southern_ocean = true,
                           south_boundary = 65.0,
                           north_boundary = 91.0,
                           barriers = nothing,
                           seed_points = [(nothing, 90.0)],
                           kw...)

    return OceanBasinMask(grid; south_boundary, north_boundary, barriers, seed_points, kw...)
end
