using Oceananigans.Grids: φnode, inactive_cell, on_architecture
using Oceananigans.Operators: Δxᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶜᶠᶜ, Δzᶠᶜᶜ
using Oceananigans.Utils: launch!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid
using KernelAbstractions: @index, @kernel

# atlantic_ocean_mask is imported by the parent module (Diagnostics.jl)

#####
##### LatitudinalBandTags struct
#####

"""
    LatitudinalBandTags{TC, G, M}

A structure for tagging grid cells with latitudinal band numbers for computing
meridional overturning circulation streamfunctions.

The tag field contains:
- `-1`: untagged (land or outside basin)
- `0, 1, 2, ...`: band number (distance from starting latitude)

This is used to identify which velocity components contribute to transport
across each latitudinal band.
"""
struct LatitudinalBandTags{TC, G, M, FT}
    tag :: TC              # Field{Center, Center, Nothing} - band number at cell centers
    grid :: G
    basin_mask :: M        # Optional OceanBasinMask or mask field
    starting_latitude :: FT
    direction :: Symbol    # :northward or :southward
end

Base.summary(lbt::LatitudinalBandTags) = "LatitudinalBandTags"

function Base.show(io::IO, lbt::LatitudinalBandTags)
    print(io, summary(lbt), '\n')
    print(io, "├── grid: ", summary(lbt.grid), '\n')
    print(io, "├── starting_latitude: ", lbt.starting_latitude, "°", '\n')
    print(io, "├── direction: ", lbt.direction, '\n')
    print(io, "└── basin_mask: ", isnothing(lbt.basin_mask) ? "nothing" : summary(lbt.basin_mask))
end

# Forward getindex to tag field
Base.getindex(lbt::LatitudinalBandTags, i, j, k) = lbt.tag[i, j, k]

#####
##### Constructor
#####

"""
    LatitudinalBandTags(grid;
                        basin_mask = nothing,
                        starting_latitude = 90.0,
                        direction = :southward)

Create latitudinal band tags for computing meridional overturning circulation.

The algorithm uses wavefront propagation starting from `starting_latitude` and
propagating in the specified `direction` (:northward or :southward). Each cell
is tagged with a band number representing its distance from the starting latitude.

Arguments
=========
- `grid`: The ocean grid (LatitudeLongitudeGrid, TripolarGrid, or ImmersedBoundaryGrid)

Keyword Arguments
=================
- `basin_mask`: Optional mask restricting tagging to a specific ocean basin.
                Can be an OceanBasinMask or any field-like object where
                `basin_mask[i, j, 1] > 0` indicates a valid cell.
- `starting_latitude`: Latitude to start the wavefront propagation. Default: 90.0
- `direction`: Direction of propagation, either `:northward` or `:southward`. Default: :southward

Returns
=======
A `LatitudinalBandTags` with a 2D tag field containing band numbers.

Example
=======
```julia
using ClimaOcean

grid = LatitudeLongitudeGrid(size=(360, 180, 50), latitude=(-90, 90), longitude=(0, 360), z=(-5000, 0))
atlantic = atlantic_ocean_mask(grid)
tags = LatitudinalBandTags(grid; basin_mask=atlantic, starting_latitude=65.0, direction=:southward)
```
"""
function LatitudinalBandTags(grid;
                              basin_mask = nothing,
                              starting_latitude = 90.0,
                              direction = :southward)

    if direction ∉ (:northward, :southward)
        throw(ArgumentError("direction must be :northward or :southward, got $direction"))
    end

    FT = eltype(grid)
    starting_latitude = convert(FT, starting_latitude)

    # Create tag field (2D, at cell centers)
    tag = Field{Center, Center, Nothing}(grid, Int32)
    set!(tag, -1)  # Initialize all cells as untagged

    tags = LatitudinalBandTags(tag, grid, basin_mask, starting_latitude, direction)

    compute_latitudinal_tags!(tags)

    return tags
end

#####
##### Wavefront propagation algorithm
#####

"""
    compute_latitudinal_tags!(tags::LatitudinalBandTags)

Compute latitudinal band tags using wavefront propagation.

The algorithm:
1. Initialize starting band (cells near starting_latitude) with tag 0
2. Iteratively propagate: for each untagged wet cell adjacent to current band,
   tag with next band number
3. Handle periodic longitude boundaries and tripolar fold
4. Continue until no more cells can be tagged
"""
function compute_latitudinal_tags!(tags::LatitudinalBandTags)
    grid = tags.grid
    arch = architecture(grid)

    # Initialize starting band
    initialize_starting_band!(tags)
    fill_halo_regions!(tags.tag)

    # Create a temporary field for double-buffering
    tag_next = Field{Center, Center, Nothing}(grid, Int32)
    parent(tag_next) .= parent(tags.tag)

    # Iterative wavefront propagation
    current_band = 0
    max_iterations = sum(size(grid)[1:2])  # Upper bound on iterations
    cells_tagged = 1  # Start with nonzero to enter loop

    while cells_tagged > 0 && current_band < max_iterations
        # Propagate to next band
        launch!(arch, grid, :xy, _propagate_wavefront!,
                tag_next, tags.tag, grid, tags.basin_mask,
                Int32(current_band), tags.direction)

        # Handle tripolar fold if applicable
        propagate_across_tripolar_fold!(tag_next, tags.tag, grid, Int32(current_band))

        fill_halo_regions!(tag_next)

        # Count newly tagged cells
        cells_tagged = count_newly_tagged(tag_next, tags.tag, grid)

        # Swap buffers
        parent(tags.tag) .= parent(tag_next)

        current_band += 1
    end

    fill_halo_regions!(tags.tag)

    return tags
end

"""
    initialize_starting_band!(tags::LatitudinalBandTags)

Initialize cells at the starting latitude with band number 0.
"""
function initialize_starting_band!(tags::LatitudinalBandTags)
    grid = tags.grid
    arch = architecture(grid)

    # Compute latitude spacing for tolerance
    Ny = size(grid, 2)
    FT = eltype(grid)
    Δφ = convert(FT, 180) / Ny

    # Convert direction to boolean for kernel
    is_southward = tags.direction == :southward

    launch!(arch, grid, :xy, _initialize_starting_band!,
            tags.tag, grid, tags.starting_latitude, Δφ,
            tags.basin_mask, is_southward)

    return nothing
end

@kernel function _initialize_starting_band!(tag, grid, starting_latitude, Δφ, basin_mask, is_southward)
    i, j = @index(Global, NTuple)

    φ = φnode(i, j, 1, grid, Center(), Center(), Center())

    # Check if cell is in basin (if mask provided)
    in_basin = isnothing(basin_mask) || basin_mask[i, j, 1] > 0

    # Check if cell is wet (not land)
    is_wet = !inactive_cell(i, j, 1, grid)

    # Check if cell is at starting latitude
    # For southward: φ >= starting_latitude - Δφ
    # For northward: φ <= starting_latitude + Δφ
    is_starting_south = φ >= starting_latitude - Δφ
    is_starting_north = φ <= starting_latitude + Δφ
    is_starting = ifelse(is_southward, is_starting_south, is_starting_north)

    @inbounds tag[i, j, 1] = ifelse(in_basin & is_wet & is_starting, Int32(0), Int32(-1))
end

@kernel function _propagate_wavefront!(tag_next, tag, grid, basin_mask, current_band, direction)
    i, j = @index(Global, NTuple)

    Nx, Ny = size(grid, 1), size(grid, 2)

    # Get current cell's tag
    @inbounds current_tag = tag[i, j, 1]

    # Check if cell is in basin (if mask provided)
    in_basin = isnothing(basin_mask) || basin_mask[i, j, 1] > 0

    # Check if cell is wet
    is_wet = !inactive_cell(i, j, 1, grid)

    # Check neighbors for current_band tag
    # Periodic in longitude (i direction)
    im1 = ifelse(i == 1, Nx, i - 1)
    ip1 = ifelse(i == Nx, 1, i + 1)

    # Not periodic in latitude (j direction) for standard grids
    jm1 = max(1, j - 1)
    jp1 = min(Ny, j + 1)

    @inbounds begin
        tag_west  = tag[im1, j, 1]
        tag_east  = tag[ip1, j, 1]
        tag_south = tag[i, jm1, 1]
        tag_north = tag[i, jp1, 1]
    end

    # Check if any neighbor has current_band tag
    neighbor_has_current = (tag_west == current_band) |
                           (tag_east == current_band) |
                           (tag_south == current_band) |
                           (tag_north == current_band)

    # Determine the new tag value:
    # - If already tagged (>= 0), keep the tag
    # - If not in basin or not wet, set to -1
    # - If neighbor has current_band, set to current_band + 1
    # - Otherwise, set to -1
    already_tagged = current_tag >= 0
    valid_cell = in_basin & is_wet
    next_band = current_band + Int32(1)

    new_tag = ifelse(already_tagged, current_tag,
              ifelse(valid_cell & neighbor_has_current, next_band, Int32(-1)))

    @inbounds tag_next[i, j, 1] = new_tag
end

#####
##### Tripolar grid handling
#####

# Default: no tripolar fold handling
propagate_across_tripolar_fold!(tag_next, tag, grid, current_band) = nothing

# Specialization for TripolarGrid
function propagate_across_tripolar_fold!(tag_next, tag, grid::TripolarGridOfSomeKind, current_band)
    arch = architecture(grid)
    launch!(arch, grid, :x, _propagate_tripolar_fold!, tag_next, tag, grid, current_band)
    return nothing
end

@kernel function _propagate_tripolar_fold!(tag_next, tag, grid, current_band)
    i = @index(Global)

    Nx, Ny = size(grid, 1), size(grid, 2)

    # At j = Ny, cells connect across the fold
    # Cell (i, Ny) connects to cell (Nx - i + 1, Ny)
    i_fold = Nx - i + 1

    @inbounds begin
        tag_here = tag[i, Ny, 1]
        tag_fold = tag[i_fold, Ny, 1]
        current_tag_next = tag_next[i, Ny, 1]
    end

    # If this cell is untagged but the fold partner has current_band
    should_tag = (tag_here < 0) & (tag_fold == current_band)
    new_tag = ifelse(should_tag, current_band + Int32(1), current_tag_next)

    @inbounds tag_next[i, Ny, 1] = new_tag
end

#####
##### Utility functions
#####

"""
    count_newly_tagged(tag_next, tag, grid)

Count cells that were newly tagged in this iteration.
"""
function count_newly_tagged(tag_next, tag, grid)
    # Count cells that changed from -1 to >= 0
    Nx, Ny = size(grid, 1), size(grid, 2)
    count = 0

    # Use field indexing which accounts for halos
    for j in 1:Ny
        for i in 1:Nx
            if tag[i, j, 1] < 0 && tag_next[i, j, 1] >= 0
                count += 1
            end
        end
    end

    return count
end

#####
##### Velocity band-crossing functions
#####

"""
    v_crosses_band(tags::LatitudinalBandTags, i, j, band)

Check if v-velocity at (i, j) crosses the boundary of latitudinal `band`.

V-velocity at (i, j) is located at (Center, Face) and represents transport
between cells (i, j-1) and (i, j). It crosses a band boundary if these
cells have different band tags.
"""
function v_crosses_band(tags::LatitudinalBandTags, i, j, band)
    @inbounds begin
        tag_north = tags.tag[i, j, 1]
        tag_south = tags.tag[i, j-1, 1]
    end

    # V crosses band if it connects cells in adjacent bands
    crosses_northward = (tag_north == band) && (tag_south == band - 1)
    crosses_southward = (tag_south == band) && (tag_north == band - 1)

    return crosses_northward || crosses_southward
end

"""
    u_crosses_band(tags::LatitudinalBandTags, i, j, band)

Check if u-velocity at (i, j) crosses the boundary of latitudinal `band`.

U-velocity at (i, j) is located at (Face, Center) and represents transport
between cells (i-1, j) and (i, j). It crosses a band boundary if these
cells have different band tags.
"""
function u_crosses_band(tags::LatitudinalBandTags, i, j, band)
    @inbounds begin
        tag_east = tags.tag[i, j, 1]
        tag_west = tags.tag[i-1, j, 1]
    end

    # U crosses band if it connects cells in adjacent bands
    crosses_eastward = (tag_east == band) && (tag_west == band - 1)
    crosses_westward = (tag_west == band) && (tag_east == band - 1)

    return crosses_eastward || crosses_westward
end

"""
    v_band_sign(tags::LatitudinalBandTags, i, j, band)

Return the sign of v-velocity contribution at (i, j) to band transport.
Returns +1 if velocity points from lower to higher band number,
-1 if from higher to lower, and 0 if not crossing the band.
"""
function v_band_sign(tags::LatitudinalBandTags, i, j, band)
    @inbounds begin
        tag_north = tags.tag[i, j, 1]
        tag_south = tags.tag[i, j-1, 1]
    end

    if (tag_north == band) && (tag_south == band - 1)
        return 1   # Positive v brings water into higher band
    elseif (tag_south == band) && (tag_north == band - 1)
        return -1  # Positive v takes water out of higher band
    else
        return 0   # Does not cross this band
    end
end

"""
    u_band_sign(tags::LatitudinalBandTags, i, j, band)

Return the sign of u-velocity contribution at (i, j) to band transport.
Returns +1 if velocity points from lower to higher band number,
-1 if from higher to lower, and 0 if not crossing the band.
"""
function u_band_sign(tags::LatitudinalBandTags, i, j, band)
    @inbounds begin
        tag_east = tags.tag[i, j, 1]
        tag_west = tags.tag[i-1, j, 1]
    end

    if (tag_east == band) && (tag_west == band - 1)
        return 1   # Positive u brings water into higher band
    elseif (tag_west == band) && (tag_east == band - 1)
        return -1  # Positive u takes water out of higher band
    else
        return 0   # Does not cross this band
    end
end

#####
##### Band transport computation
#####

"""
    compute_band_transport(u, v, tags::LatitudinalBandTags, band_number)

Compute the meridional volume transport across a specific latitudinal band.

Returns a 1D array of transport at each depth level (positive = increasing band number).
"""
function compute_band_transport(u, v, tags::LatitudinalBandTags, band_number)
    grid = tags.grid
    arch = architecture(grid)
    Nz = size(grid, 3)

    # Transport at each depth level
    transport = zeros(eltype(grid), Nz)
    transport = on_architecture(arch, transport)

    # GPU kernel would need atomic operations or reduction
    # For now, implement CPU version
    _compute_band_transport_cpu!(transport, u, v, tags, grid, Int32(band_number))

    return Array(transport)
end

function _compute_band_transport_cpu!(transport, u, v, tags, grid, band_number)
    Nx, Ny, Nz = size(grid)

    for k in 1:Nz
        total = zero(eltype(grid))

        # V-velocity contributions
        for j in 2:Ny  # j=1 has no j-1
            for i in 1:Nx
                sign = v_band_sign(tags, i, j, band_number)
                if sign != 0
                    # Transport = v * Δx * Δz
                    Δx = Δxᶜᶠᶜ(i, j, k, grid)
                    Δz = Δzᶜᶠᶜ(i, j, k, grid)
                    total += sign * v[i, j, k] * Δx * Δz
                end
            end
        end

        # U-velocity contributions
        for j in 1:Ny
            for i in 2:Nx  # i=1 needs periodic handling
                sign = u_band_sign(tags, i, j, band_number)
                if sign != 0
                    # Transport = u * Δy * Δz
                    Δy = Δyᶠᶜᶜ(i, j, k, grid)
                    Δz = Δzᶠᶜᶜ(i, j, k, grid)
                    total += sign * u[i, j, k] * Δy * Δz
                end
            end
            # Handle i=1 (periodic)
            sign = u_band_sign(tags, 1, j, band_number)
            if sign != 0
                Δy = Δyᶠᶜᶜ(1, j, k, grid)
                Δz = Δzᶠᶜᶜ(1, j, k, grid)
                total += sign * u[1, j, k] * Δy * Δz
            end
        end

        transport[k] = total
    end

    return nothing
end

#####
##### AMOC streamfunction
#####

"""
    compute_amoc_streamfunction(u, v, tags::LatitudinalBandTags)

Compute the Atlantic Meridional Overturning Circulation streamfunction.

Returns a 2D array (n_bands × Nz) representing ψ(φ, z), where:
- ψ[band, k] represents the cumulative transport from the surface to depth level k
- Positive values indicate clockwise overturning (in standard view)

The streamfunction is computed by:
1. Computing transport across each band boundary at each depth
2. Cumulatively summing from the surface downward
"""
function compute_amoc_streamfunction(u, v, tags::LatitudinalBandTags)
    grid = tags.grid
    Nx, Ny, Nz = size(grid)

    # Find max band number using field indexing
    max_band = 0
    for j in 1:Ny
        for i in 1:Nx
            band = tags.tag[i, j, 1]
            max_band = max(max_band, band)
        end
    end

    if max_band < 1
        error("No valid latitudinal bands found. Check basin_mask and starting_latitude.")
    end

    # Initialize streamfunction array
    FT = eltype(grid)
    ψ = zeros(FT, max_band, Nz)

    # Compute transport for each band
    for band in 1:max_band
        transport = compute_band_transport(u, v, tags, band)

        # Cumulative sum from surface (k=Nz) downward
        ψ[band, Nz] = transport[Nz]
        for k in (Nz-1):-1:1
            ψ[band, k] = ψ[band, k+1] + transport[k]
        end
    end

    return ψ
end

#####
##### Band latitude mapping
#####

"""
    band_latitudes(tags::LatitudinalBandTags)

Return the approximate latitude corresponding to each band number.

Computes the mean latitude of all cells in each band.
"""
function band_latitudes(tags::LatitudinalBandTags)
    grid = tags.grid
    Nx, Ny = size(grid, 1), size(grid, 2)

    # Find max band using field indexing
    max_band = 0
    for j in 1:Ny
        for i in 1:Nx
            band = tags.tag[i, j, 1]
            max_band = max(max_band, band)
        end
    end

    if max_band < 1
        return Float64[]
    end

    # Sum latitudes and count cells for each band
    lat_sum = zeros(max_band)
    lat_count = zeros(Int, max_band)

    for j in 1:Ny
        for i in 1:Nx
            band = tags.tag[i, j, 1]
            if band >= 1  # Skip untagged (-1) and starting band (0)
                φ = φnode(i, j, 1, grid, Center(), Center(), Center())
                lat_sum[band] += φ
                lat_count[band] += 1
            end
        end
    end

    # Compute mean latitudes
    latitudes = similar(lat_sum)
    for band in 1:max_band
        latitudes[band] = lat_count[band] > 0 ? lat_sum[band] / lat_count[band] : NaN
    end

    return latitudes
end

#####
##### Convenience functions
#####

"""
    atlantic_latitudinal_bands(grid;
                                south_boundary = -34.0,
                                north_boundary = 65.0)

Create latitudinal band tags for the Atlantic Ocean basin.

This is a convenience function that creates an Atlantic basin mask and
sets up appropriate tagging for AMOC computation.

Arguments
=========
- `grid`: The ocean grid

Keyword Arguments
=================
- `south_boundary`: Southern boundary of Atlantic mask. Default: -34.0
- `north_boundary`: Northern boundary of Atlantic mask. Default: 65.0

Returns
=======
A `LatitudinalBandTags` configured for Atlantic AMOC computation.
"""
function atlantic_latitudinal_bands(grid;
                                     south_boundary = -34.0,
                                     north_boundary = 65.0)
    atlantic = atlantic_ocean_mask(grid; south_boundary, north_boundary)
    return LatitudinalBandTags(grid;
                                basin_mask = atlantic,
                                starting_latitude = north_boundary,
                                direction = :southward)
end
