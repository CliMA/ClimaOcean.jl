using Oceananigans.Grids: φnode, λnode, inactive_cell, on_architecture, znode
using Oceananigans.Operators: Δxᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶜᶠᶜ, Δzᶠᶜᶜ
using Oceananigans.Fields: interior
using Oceananigans.Utils: launch!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid, TripolarGridOfSomeKind
using KernelAbstractions: @index, @kernel

# atlantic_ocean_mask is imported by the parent module (Diagnostics.jl)

#####
##### IsoLatitudeBrokenLine - MITgcm-style broken line representation
#####

"""
    IsoLatitudeBrokenLine{I, F}

Represents a broken line approximating an iso-latitude contour on a curvilinear grid.

On non-lat-lon grids (tripolar, cubed-sphere), latitude circles don't align with
grid cell edges. The broken line is constructed by finding all velocity points
(U and V) that separate cells in different latitude bands.

Fields
======
- `latitude`: Target latitude for this broken line
- `n_points`: Number of velocity points on this broken line
- `i_indices`: i-indices of velocity points (Vector{Int})
- `j_indices`: j-indices of velocity points (Vector{Int})
- `u_factor`: Weight for U-velocity: +1, -1, or 0 (Vector{Int8})
- `v_factor`: Weight for V-velocity: +1, -1, or 0 (Vector{Int8})

The factors encode both component and sign:
- `u_factor = +1`: U-point, positive U = northward transport
- `u_factor = -1`: U-point, negative U = northward transport
- `v_factor = +1`: V-point, positive V = northward transport
- `v_factor = -1`: V-point, negative V = northward transport
- Factor = 0 means this velocity component doesn't contribute at this point
"""
struct IsoLatitudeBrokenLine{FT, I, F}
    latitude  :: FT
    n_points  :: Int
    i_indices :: I
    j_indices :: I
    u_factor  :: F
    v_factor  :: F
end

Base.summary(bl::IsoLatitudeBrokenLine) = "IsoLatitudeBrokenLine at $(bl.latitude)° with $(bl.n_points) points"

#####
##### BrokenLineSet - collection of broken lines for multiple latitudes
#####

"""
    BrokenLineSet{FT, BL, TC, G, M}

A collection of iso-latitude broken lines for computing meridional overturning
streamfunctions using the MITgcm algorithm.

This structure precomputes the broken line geometry once from the grid, then
reuses it for all transport calculations.

Fields
======
- `latitudes`: Vector of target latitudes
- `lines`: Vector of IsoLatitudeBrokenLine objects
- `tag`: 2D field of latitude band tags at cell centers
- `grid`: The underlying grid
- `basin_mask`: Optional mask for restricting to a specific basin
"""
struct BrokenLineSet{FT, BL, TC, G, M}
    latitudes  :: Vector{FT}
    lines      :: BL
    tag        :: TC
    grid       :: G
    basin_mask :: M
end

Base.summary(bls::BrokenLineSet) = "BrokenLineSet with $(length(bls.latitudes)) latitudes"

function Base.show(io::IO, bls::BrokenLineSet)
    print(io, summary(bls), '\n')
    print(io, "├── grid: ", summary(bls.grid), '\n')
    print(io, "├── latitude range: ", minimum(bls.latitudes), "° to ", maximum(bls.latitudes), "°", '\n')
    print(io, "└── basin_mask: ", isnothing(bls.basin_mask) ? "nothing" : summary(bls.basin_mask))
end

#####
##### Constructor
#####

"""
    BrokenLineSet(grid, latitudes;
                  basin_mask = nothing)

Create a set of iso-latitude broken lines for computing meridional overturning circulation.

This implements the MITgcm algorithm from `mk_isoLat_bkl.m`:
1. Tag each cell center with a latitude band number based on its actual latitude
2. For each target latitude, find all velocity points (U and V) that separate
   cells in different latitude bands
3. Store the sign indicating whether positive velocity = northward transport

Arguments
=========
- `grid`: The ocean grid (LatitudeLongitudeGrid, TripolarGrid, or ImmersedBoundaryGrid)
- `latitudes`: Vector of target latitudes for broken lines

Keyword Arguments
=================
- `basin_mask`: Optional mask restricting computation to a specific ocean basin.
                Can be an OceanBasinMask or any field-like object where
                `basin_mask[i, j, 1] > 0` indicates a valid cell.

Returns
=======
A `BrokenLineSet` containing precomputed broken line geometry.

Example
=======
```julia
using ClimaOcean

grid = TripolarGrid(size=(360, 180, 50), z=(-5000, 0))
atlantic = atlantic_ocean_mask(grid)
latitudes = collect(-34.0:1.0:65.0)
broken_lines = BrokenLineSet(grid, latitudes; basin_mask=atlantic)
ψ = compute_streamfunction(u, v, broken_lines)
```
"""
function BrokenLineSet(grid, latitudes;
                       basin_mask = nothing)

    FT = eltype(grid)
    latitudes = convert.(FT, latitudes)

    # Sort latitudes for consistent ordering
    latitudes = sort(latitudes)

    # Step 1: Tag cell centers by latitude band
    tag = compute_latitude_tags(grid, latitudes, basin_mask)

    # Step 2: Generate broken lines for each latitude
    lines = generate_broken_lines(grid, latitudes, tag)

    return BrokenLineSet(latitudes, lines, tag, grid, basin_mask)
end

#####
##### Latitude band tagging (MITgcm Step 1)
#####

"""
    compute_latitude_tags(grid, latitudes, basin_mask)

Tag each cell center with a latitude band number.

Cell (i,j) gets tag `jl` if `latitudes[jl] <= φ[i,j] < latitudes[jl+1]`.
Land cells and cells outside the basin mask get tag -1.
Cells south of all latitudes get tag 0.
Cells north of all latitudes get tag `length(latitudes)`.
"""
function compute_latitude_tags(grid, latitudes, basin_mask)
    Nx, Ny, Nz = size(grid)
    FT = eltype(grid)
    n_lats = length(latitudes)

    # Create tag field (2D, at cell centers)
    tag = Field{Center, Center, Nothing}(grid)

    arch = architecture(grid)

    launch!(arch, grid, :xy, _compute_latitude_tags!,
            tag, grid, latitudes, n_lats, basin_mask)

    fill_halo_regions!(tag)

    return tag
end

@kernel function _compute_latitude_tags!(tag, grid, latitudes, n_lats, basin_mask)
    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)

    # Get cell center latitude
    φ = φnode(i, j, 1, grid, Center(), Center(), Center())

    # Check if cell is in basin (if mask provided)
    in_basin = isnothing(basin_mask) || basin_mask[i, j, 1] > 0

    # Check if cell is wet (not land) - use surface level (Nz) since
    # deeper levels may be underground in shallow regions
    is_wet = !inactive_cell(i, j, Nz, grid)

    # Check for valid latitude (guard against NaN at poles or special grid points)
    valid_φ = isfinite(φ)

    if !in_basin || !is_wet || !valid_φ
        @inbounds tag[i, j, 1] = -1
    else
        # Find which latitude band this cell belongs to
        # Tag = jl means latitudes[jl] <= φ < latitudes[jl+1]
        # Tag = 0 means φ < latitudes[1]
        # Tag = n_lats means φ >= latitudes[n_lats]

        cell_tag = 0
        for jl in 1:n_lats
            if φ >= latitudes[jl]
                cell_tag = jl
            end
        end

        @inbounds tag[i, j, 1] = cell_tag
    end
end

#####
##### Broken line generation (MITgcm Step 2)
#####

"""
    generate_broken_lines(grid, latitudes, tag)

Generate broken lines for each target latitude.

A velocity point lies on the broken line for latitude `jl` if it separates
two cells whose tags straddle `jl`:
- For U-points between (i-1,j) and (i,j): crosses if tags straddle jl
- For V-points between (i,j-1) and (i,j): crosses if tags straddle jl

The sign indicates whether positive velocity = northward transport.
"""
function generate_broken_lines(grid, latitudes, tag)
    Nx, Ny, Nz = size(grid)
    n_lats = length(latitudes)

    # Move tag to CPU for line generation and get the interior data
    tag_cpu = on_architecture(CPU(), tag)
    # Use interior to ensure we get the actual data array
    tag_data = interior(tag_cpu, :, :, 1)

    lines = IsoLatitudeBrokenLine[]

    for (jl, lat) in enumerate(latitudes)
        i_indices = Int[]
        j_indices = Int[]
        u_factors = Int8[]
        v_factors = Int8[]

        # Check all U-points (between cells i-1,j and i,j)
        # U-point at (i,j) is at Face location in x, Center in y
        for j in 1:Ny
            for i in 1:Nx
                # Handle periodic boundary: i-1 wraps to Nx
                im1 = i == 1 ? Nx : i - 1

                tag_east = tag_data[i, j]
                tag_west = tag_data[im1, j]

                # Skip if either cell is land (tag = -1)
                # But allow cells outside the latitude range (tag = 0 or tag = n_lats)
                if tag_east == -1 || tag_west == -1
                    continue
                end

                # U-point crosses latitude jl if cells straddle jl
                # tag >= jl means cell is at or north of latitude jl
                # tag < jl means cell is south of latitude jl

                if tag_east >= jl && tag_west < jl
                    # Positive U (eastward) moves water from west to east
                    # West cell is south of jl, east cell is at/north of jl
                    # So positive U contributes to northward transport
                    push!(i_indices, i)
                    push!(j_indices, j)
                    push!(u_factors, Int8(+1))
                    push!(v_factors, Int8(0))
                elseif tag_east < jl && tag_west >= jl
                    # Positive U moves water from at/north of jl to south of jl
                    # So positive U contributes to southward transport (negative sign)
                    push!(i_indices, i)
                    push!(j_indices, j)
                    push!(u_factors, Int8(-1))
                    push!(v_factors, Int8(0))
                end
            end
        end

        # Check all V-points (between cells i,j-1 and i,j)
        # V-point at (i,j) is at Center in x, Face in y
        for j in 2:Ny  # j=1 has no j-1 neighbor (southern boundary)
            for i in 1:Nx
                tag_north = tag_data[i, j]
                tag_south = tag_data[i, j-1]

                # Skip if either cell is land (tag = -1)
                if tag_north == -1 || tag_south == -1
                    continue
                end

                # V-point crosses latitude jl if cells straddle jl
                if tag_north >= jl && tag_south < jl
                    # Positive V (northward) moves water from south to north
                    # South cell is below jl, north cell is at/above jl
                    # So positive V = northward transport
                    push!(i_indices, i)
                    push!(j_indices, j)
                    push!(u_factors, Int8(0))
                    push!(v_factors, Int8(+1))
                elseif tag_north < jl && tag_south >= jl
                    # Positive V moves water from at/above jl to below jl
                    # So positive V = southward transport (negative sign)
                    push!(i_indices, i)
                    push!(j_indices, j)
                    push!(u_factors, Int8(0))
                    push!(v_factors, Int8(-1))
                end
            end
        end

        # Handle tripolar fold at j = Ny
        add_tripolar_fold_points!(i_indices, j_indices, u_factors, v_factors,
                                  grid, tag_data, jl, Ny)

        n_points = length(i_indices)

        push!(lines, IsoLatitudeBrokenLine(lat, n_points,
                                           i_indices, j_indices,
                                           u_factors, v_factors))
    end

    return lines
end

# Default: no tripolar fold handling
add_tripolar_fold_points!(i_indices, j_indices, u_factors, v_factors, grid, tag_data, jl, Ny) = nothing

# Specialization for TripolarGrid
function add_tripolar_fold_points!(i_indices, j_indices, u_factors, v_factors,
                                   grid::TripolarGridOfSomeKind, tag_data, jl, Ny)
    Nx = size(grid, 1)

    # At j = Ny, cells connect across the fold
    # Cell (i, Ny) connects to cell (Nx - i + 1, Ny) across the fold
    # This creates additional V-like connections at the fold

    for i in 1:Nx÷2  # Only check half to avoid double-counting
        i_fold = Nx - i + 1

        tag_here = tag_data[i, Ny]
        tag_fold = tag_data[i_fold, Ny]

        # Skip if either cell is land
        if tag_here == -1 || tag_fold == -1
            continue
        end

        # Check if the fold connection crosses latitude jl
        if tag_here >= jl && tag_fold < jl
            push!(i_indices, i)
            push!(j_indices, Ny)
            push!(u_factors, Int8(0))
            push!(v_factors, Int8(+1))  # Convention: positive = toward higher tag
        elseif tag_here < jl && tag_fold >= jl
            push!(i_indices, i)
            push!(j_indices, Ny)
            push!(u_factors, Int8(0))
            push!(v_factors, Int8(-1))
        end
    end

    return nothing
end

#####
##### Angle computation for grid-to-geographic transformation
#####

"""
    compute_face_meridional_angle(grid, i, j, is_u_face)

Compute the cosine of the angle between a face normal and the northward direction.

For curvilinear grids (tripolar, cubed-sphere), the grid axes are not aligned with
geographic directions. To compute true meridional transport, we need to project
the face-normal transport onto the northward direction.

Returns cos(θ) where θ is the angle between the face normal and north.
- cos(θ) ≈ 1 for V-faces on a lat-lon grid (face normal is northward)
- cos(θ) ≈ 0 for U-faces on a lat-lon grid (face normal is eastward)

For `is_u_face=true`:  U-face at (i, j), normal points from cell (i-1, j) to (i, j)
For `is_u_face=false`: V-face at (i, j), normal points from cell (i, j-1) to (i, j)
"""
function compute_face_meridional_angle(grid, i, j, is_u_face)
    Nx, Ny = size(grid, 1), size(grid, 2)

    if is_u_face
        # U-face: normal from (i-1, j) to (i, j)
        im1 = i == 1 ? Nx : i - 1
        λ1 = λnode(im1, j, 1, grid, Center(), Center(), Center())
        φ1 = φnode(im1, j, 1, grid, Center(), Center(), Center())
        λ2 = λnode(i, j, 1, grid, Center(), Center(), Center())
        φ2 = φnode(i, j, 1, grid, Center(), Center(), Center())
    else
        # V-face: normal from (i, j-1) to (i, j)
        jm1 = max(1, j - 1)
        λ1 = λnode(i, jm1, 1, grid, Center(), Center(), Center())
        φ1 = φnode(i, jm1, 1, grid, Center(), Center(), Center())
        λ2 = λnode(i, j, 1, grid, Center(), Center(), Center())
        φ2 = φnode(i, j, 1, grid, Center(), Center(), Center())
    end

    # Handle NaN coordinates
    if !isfinite(λ1) || !isfinite(φ1) || !isfinite(λ2) || !isfinite(φ2)
        return is_u_face ? 0.0 : 1.0  # Default: U-faces are zonal, V-faces are meridional
    end

    # Compute direction in local Cartesian coordinates
    # East component: proportional to Δλ * cos(φ_avg)
    # North component: proportional to Δφ
    Δλ = λ2 - λ1
    Δφ = φ2 - φ1
    φ_avg = 0.5 * (φ1 + φ2)

    # Handle periodic longitude (Δλ should be small for adjacent cells)
    if Δλ > 180
        Δλ -= 360
    elseif Δλ < -180
        Δλ += 360
    end

    # Convert to approximate Cartesian (east, north) components
    # Using cos(φ) to account for converging meridians
    Δeast = Δλ * cosd(φ_avg)
    Δnorth = Δφ

    # Compute magnitude and angle
    magnitude = sqrt(Δeast^2 + Δnorth^2)

    if magnitude < 1e-10
        return is_u_face ? 0.0 : 1.0  # Default for degenerate cases
    end

    # cos(angle from north) = Δnorth / magnitude
    return Δnorth / magnitude
end

#####
##### Transport computation (MITgcm Step 4 with angle correction)
#####

"""
    compute_volume_transport(u, v, grid, broken_line, k)

Compute meridional volume transport across a broken line at depth level k.

This includes proper angle correction for curvilinear grids (tripolar, cubed-sphere).
On these grids, the velocity components (u, v) are in grid coordinates, not
geographic coordinates. To get true meridional transport:

- For V-faces: transport_meridional = v * Δx * Δz * cos(θ_v)
- For U-faces: transport_meridional = u * Δy * Δz * cos(θ_u)

where θ is the angle between the face normal and the northward direction.

Returns the total meridional volume transport in m³/s (positive = northward).
"""
function compute_volume_transport(u, v, grid, broken_line::IsoLatitudeBrokenLine, k)
    transport = 0.0  # Use Float64 explicitly
    Nx = size(grid, 1)

    for n in 1:broken_line.n_points
        i = broken_line.i_indices[n]
        j = broken_line.j_indices[n]
        uf = broken_line.u_factor[n]
        vf = broken_line.v_factor[n]

        if uf != 0
            # U-face transport with angle correction
            u_val = u[i, j, k]
            Δy = Δyᶠᶜᶜ(i, j, k, grid)
            Δz = Δzᶠᶜᶜ(i, j, k, grid)

            if isfinite(u_val) && isfinite(Δy) && isfinite(Δz) && Δz > 0
                # Angle correction: project onto northward direction
                cos_angle = compute_face_meridional_angle(grid, i, j, true)
                transport += uf * u_val * Δy * Δz * cos_angle
            end
        end

        if vf != 0
            # V-face transport with angle correction
            v_val = v[i, j, k]
            Δx = Δxᶜᶠᶜ(i, j, k, grid)
            Δz = Δzᶜᶠᶜ(i, j, k, grid)

            if isfinite(v_val) && isfinite(Δx) && isfinite(Δz) && Δz > 0
                # Angle correction: project onto northward direction
                cos_angle = compute_face_meridional_angle(grid, i, j, false)
                transport += vf * v_val * Δx * Δz * cos_angle
            end
        end
    end

    return transport
end

#####
##### Streamfunction computation (MITgcm Step 5)
#####

"""
    compute_streamfunction(u, v, broken_lines::BrokenLineSet)

Compute the meridional overturning streamfunction Ψ(φ, z).

The streamfunction is defined at vertical faces (Nz+1 values) and computed by
integrating transport from bottom to surface:

    Ψ(z) = ∫_{-H}^{z} v̄ dz'

Integration proceeds from bottom (k=1) upward to surface (k=Nz):

    Ψ[k+1] = Ψ[k] + transport[k]

Returns a 2D array (n_latitudes × (Nz+1)) where:
- Ψ[jl, k] is the streamfunction at latitude `latitudes[jl]` and face k
- Ψ[:, 1] = 0 at the bottom (boundary condition)
- Ψ[:, Nz+1] should be ≈ 0 at the surface for mass conservation
- Units: m³/s (divide by 1e6 for Sverdrups)

For AMOC, the maximum |Ψ| typically occurs around 1000m depth with magnitude ~15-20 Sv.
"""
function compute_streamfunction(u, v, broken_lines::BrokenLineSet)
    grid = broken_lines.grid
    Nz = size(grid, 3)
    n_lats = length(broken_lines.latitudes)

    # Use Float64 explicitly to avoid potential issues with grid element types
    FT = Float64

    # Streamfunction at faces: Nz+1 values
    # ψ[jl, k] is the streamfunction at face k for latitude jl
    # Face 1 is at the bottom, face Nz+1 is at the surface
    ψ = zeros(FT, n_lats, Nz + 1)

    for (jl, line) in enumerate(broken_lines.lines)
        # Skip empty broken lines
        if line.n_points == 0
            continue
        end

        # Bottom boundary condition: ψ = 0 at the bottom face
        ψ[jl, 1] = zero(FT)

        # Integrate from bottom (k=1) up to surface (k=Nz)
        # ψ(z) = ∫_{-H}^{z} v̄ dz'
        # In discrete form: ψ[k+1] = ψ[k] + transport[k]
        for k in 1:Nz
            transport_k = compute_volume_transport(u, v, grid, line, k)
            # Guard against NaN propagation
            if !isfinite(transport_k)
                transport_k = zero(FT)
            end
            ψ[jl, k+1] = ψ[jl, k] + transport_k
        end
    end

    return ψ
end
