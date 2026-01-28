using Oceananigans.Grids: φnode, λnode, inactive_cell, on_architecture, znode
using Oceananigans.Operators: Δxᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶜᶠᶜ, Δzᶠᶜᶜ, Δzᶜᶜᶜ, extrinsic_vector
using Oceananigans.Fields: interior
using Oceananigans.Utils: launch!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid, TripolarGridOfSomeKind
using KernelAbstractions: @index, @kernel

#####
##### Latitude band tagging (MITgcm Step 1)
#####

function compute_latitude_tags(grid, latitudes)
    Nx, Ny, Nz = size(grid)
    FT = eltype(grid)
    Nlats = length(latitudes)

    # Create tag field (2D, at cell centers)
    tag  = Field{Center, Center, Nothing}(grid)

    launch!(architecture(grid), grid, :xy, _compute_latitude_tags!, 
            tag, grid, latitudes, Nlats)

    fill_halo_regions!(tag)

    return tag
end

@kernel function _compute_latitude_tags!(tag, grid, latitudes, Nlats)
    i, j = @index(Global, NTuple)
    Nz   = size(grid, 3)

    φ = φnode(i, j, 1, grid, Center(), Center(), Center())
    immersed = inactive_cell(i, j, Nz, grid)

    # Find which latitude band this cell belongs to
    # Tag = jl means latitudes[jl] <= φ < latitudes[jl+1]
    # Tag = 0 means φ < latitudes[1]
    # Tag = Nlats means φ >= latitudes[Nlats]
    cell_tag = 0
    for jl in 1:Nlats
        if φ >= latitudes[jl]
            cell_tag = jl
        end
    end

    @inbounds tag[i, j, 1] = ifelse(immersed, -1, cell_tag)
end

#####
##### Transport computation (MITgcm Step 4 with angle correction)
#####

@inline function compute_volume_transport(i, j, k, grid, u, v)
    up = extrinsic_vector(i, j, k, grid, u, ZeroField())[2]
    vp = extrinsic_vector(i, j, k, grid, ZeroField(), v)[2]
    
    Δy   = Δyᶠᶜᶜ(i, j, k, grid)
    Δx   = Δxᶜᶠᶜ(i, j, k, grid)
    ΔzFC = Δzᶠᶜᶜ(i, j, k, grid)
    ΔzCF = Δzᶜᶠᶜ(i, j, k, grid)

    return vp * Δx * ΔzCF + up * Δy * ΔzFC
end

#####
##### Streamfunction computation (MITgcm Step 5)
#####

struct TagCondition{T, M, I} <:Function
    tags :: T
    mask :: M
    target :: I
    level :: I
end

@inline (t::TagCondition)(i, j, k, grid, co) = @inbounds (t.tags[i, j, 1] == t.target) & (t.mask[i, j, 1] == 1) & (k == t.level)
@inline (t::TagCondition{<:Any, <:Nothing})(i, j, k, grid, co) = @inbounds (t.tags[i, j, 1] == t.target) & (k == t.level)

using Oceananigans.Fields: ZeroField

function compute_streamfunction(u, v, grid, latitudes; mask = nothing)
    Nz = size(grid, 3)
    Nlats = length(latitudes)
    FT = eltype(grid)

    tags = compute_latitude_tags(grid, latitudes)
    transport = KernelFunctionOperation{Face, Face, Center}(compute_volume_transport, grid, u, v)

    # Streamfunction at faces: Nz+1 values
    # ψ[jl, k] is the streamfunction at face k for latitude jl
    # Face 1 is at the bottom, face Nz+1 is at the surface
    ψ = zeros(FT, Nlats, Nz+1)

    for jl in 1:Nlats
        for k in 1:Nz
            condition = TagCondition(tags, mask, jl, k)
            ψk = sum(transport; condition)
            ψk = ifelse(!isfinite(ψk), zero(FT), ψk)
            ψ[jl, k+1] = ψ[jl, k] + ψk
            @show k, jl
        end
    end

    return ψ
end
