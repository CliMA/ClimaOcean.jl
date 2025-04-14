using ConservativeRegridding

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Grids: λnode, φnode, on_architecture, AbstractGrid
using GeometryOps: Point2

using KernelAbstractions: @kernel, @index

using SpeedyWeather.RingGrids

"""
    get_faces(grid)

Returns a list representing all horizontal grid cells in a curvilinear `grid`. 
The output is an Array of 5 * M `Point2` elements where `M = Nx * Ny`. Each column lists the vertices 
associated with a horizontal cell in clockwise order starting from the southwest (bottom left) corner.
"""
function get_faces(grid::AbstractGrid)
    Nx, Ny, _ = size(grid)
    FT = eltype(grid)

    cpu_grid = on_architecture(Oceananigans.CPU(), grid)

    sw  = fill(Point2{FT}(0, 0), 1, Nx*Ny)
    nw  = fill(Point2{FT}(0, 0), 1, Nx*Ny)
    ne  = fill(Point2{FT}(0, 0), 1, Nx*Ny)
    se  = fill(Point2{FT}(0, 0), 1, Nx*Ny)

    launch!(Oceananigans.CPU(), cpu_grid, :xy, _get_vertices!, sw, nw, ne, se, grid)
    
    vertices = vcat(sw, nw, ne, se, sw)
    
    return vertices
end

# Move to SpeedyWeather?
function get_faces(grid::SpeedyWeather.SpectralGrid)

    Grid = grid.Grid
    nlat_half = grid.nlat_half
    npoints = RingGrids.get_npoints2D(Grid, nlat_half)

    # vertex east, south, west, north (i.e. clockwise for every grid point)
    E, S, W, N = RingGrids.get_vertices(Grid, nlat_half)

    # allocate faces as Point2{Float64} so that no data copy has to be made in Makie
    faces = Matrix{NTuple{2, Float64}}(undef, 5, npoints)

    @inbounds for ij in 1:npoints
        faces[1, ij] = (E[1, ij], E[2, ij])  # clockwise
        faces[2, ij] = (S[1, ij], S[2, ij])
        faces[3, ij] = (W[1, ij], W[2, ij])
        faces[4, ij] = (N[1, ij], N[2, ij])
        faces[5, ij] = (E[1, ij], E[2, ij])  # back to east to close the polygon        
    end

    return faces
end

@kernel function _get_vertices!(sw, nw, ne, se, grid)
    i, j = @index(Global, NTuple)

    FT  = eltype(grid)
    Nx  = size(grid, 1)
    λ⁻⁻ = λnode(i,   j,   1, grid, Face(), Face(), nothing)
    λ⁺⁻ = λnode(i,   j+1, 1, grid, Face(), Face(), nothing)
    λ⁻⁺ = λnode(i+1, j,   1, grid, Face(), Face(), nothing)
    λ⁺⁺ = λnode(i+1, j+1, 1, grid, Face(), Face(), nothing)
    
    φ⁻⁻ = φnode(i,   j,   1, grid, Face(), Face(), nothing)
    φ⁺⁻ = φnode(i,   j+1, 1, grid, Face(), Face(), nothing)
    φ⁻⁺ = φnode(i+1, j,   1, grid, Face(), Face(), nothing)
    φ⁺⁺ = φnode(i+1, j+1, 1, grid, Face(), Face(), nothing)

    sw[i+(j-1)*Nx] = Point2{FT}(λ⁻⁻, φ⁻⁻)  
    nw[i+(j-1)*Nx] = Point2{FT}(λ⁻⁺, φ⁻⁺)
    ne[i+(j-1)*Nx] = Point2{FT}(λ⁺⁺, φ⁺⁺)
    se[i+(j-1)*Nx] = Point2{FT}(λ⁺⁻, φ⁺⁻)
end

