using OrthogonalSphericalShellGrids
using OrthogonalSphericalShellGrids: TRG
using Oceananigans
using Oceananigans.Grids: OSSG
using Oceananigans.Fields: fractional_index, fractional_z_index

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
@inline function fractional_horizontal_indices(λ₀, φ₀, loc, grid)
    # This is a "naive" algorithm, so it is going to be quite slow!
    # Optimizations are welcome!
    λ = λnodes(loc, grid)
    φ = φnodes(loc, grid)
    
    Nx, Ny, _ = size(grid)

    # Distance between λ₀ and a possible λ
    λd₀ = Inf

    # Initial index
    ii = 1
    jj = 1

    @inbounds begin
        for i in Nx
            # Find j-indices corresponding to φ₀
            φi = view(φ, i, :)
            j  = fractional_index(φ₀, φi, Ny) - 1
            
            j⁺ = ceil(Int, j)
            j⁻ = floor(Int, j)

            # Find the distance between λ and λ₀
            λd⁺ = norm(λ[i, j⁺] - λ₀)
            λd⁻ = norm(λ[i, j⁻] - λ₀)
            λd  = min(λd⁺, λd⁻)

            # Update the ii index depending on the distance
            # Also update previous distance
            ii  = ifelse(λd < λd₀, i, ii)
            jj  = ifelse(λd < λd₀, j, jj)
            λd₀ = ifelse(λd < λd₀, λd, λd₀)
        end

        # Check whether we want ii + 1 or ii - 1 as the second index
        λd⁺⁺ = norm(λ[ii+1, jj⁺] - λ₀)
        λd⁻⁻ = norm(λ[ii-1, jj⁻] - λ₀)
        λd⁺⁻ = norm(λ[ii+1, jj⁺] - λ₀)
        λd⁻⁺ = norm(λ[ii-1, jj⁻] - λ₀)

        ii = ifelse(λd⁺ > λd⁻, 
                    1 / (λ[ii+1, jj⁺]  - λ[ii, jj⁺] ) * (val - x₁) + i₁)
    end
   
    return ii, jj
end