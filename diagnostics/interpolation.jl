using OrthogonalSphericalShellGrids
using OrthogonalSphericalShellGrids: TRG
using Oceananigans
using Oceananigans.Grids: OSSG
using Oceananigans.Fields: fractional_z_index, index_binary_search

import Oceananigans.Fields: interpolate, fractional_indices

const ZFlatOSSG = OrthogonalSphericalShellGrid{<:Any, <:Any, <:Any, <:Flat}

# We assume OSG is at least 2-dimensional (no Flat topologies in x and y)
@inline function fractional_indices((x, y, z), grid::OSG, ℓx, ℓy, ℓz)
    ii, jj = fractional_horizontal_indices(x, y, (ℓx, ℓy, ℓz), grid)
    kk     = fractional_z_index(z, (ℓx, ℓy, ℓz), grid)
    return (ii, jj, kk)
end

# We assume OSG is at least 2-dimensional (no Flat topologies in x and y)
@inline function fractional_indices((x, y), grid::ZFlatOSSG, ℓx, ℓy, ℓz)
    ii, jj = fractional_horizontal_indices(x, y, (ℓx, ℓy, nothing), grid)
    return (ii, jj, nothing)
end

# We assume that in an OSSG, the latitude lines for a given i - index are sorted
# i.e. φ is monotone in j. This is not the case for λ that might jump between 0 and 360.
@inline function fractional_horizontal_indices(λ, φ, loc, grid)
    
end