using ConservativeRegridding

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Grids: λnode, φnode, on_architecture, AbstractGrid
using GeometryOps: Point2

using KernelAbstractions: @kernel, @index

"""
    list_cell_vertices(grid)

Returns a list representing all horizontal grid cells in a curvilinear `grid`. 
The output is an Array of 5 * M `Point2` elements where `M = Nx * Ny`. Each column lists the vertices 
associated with a horizontal cell in clockwise order starting from the southwest (bottom left) corner.
"""
function list_cell_vertices(grid)
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

