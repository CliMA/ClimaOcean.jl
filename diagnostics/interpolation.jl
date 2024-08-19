using OrthogonalSphericalShellGrids
using OrthogonalSphericalShellGrids: TRG
using Oceananigans
using Oceananigans.Grids: OSSG, λnodes, φnodes
using Oceananigans.Fields: fractional_index, fractional_z_index

import Oceananigans.Fields: interpolate, fractional_indices

ZFlatOSSG = OrthogonalSphericalShellGrid{<:Any, <:Any, <:Any, <:Flat}

# We assume OSG is at least 2-dimensional (no Flat topologies in x and y)
@inline function fractional_indices((x, y, z), grid::OSSG, ℓx, ℓy, ℓz)
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
    λ = λnodes(grid, loc...)
    φ = φnodes(grid, loc...)
    
    Nx, Ny, _ = size(grid)

    φi = view(φ, 1, :)
    j₁ = fractional_index(φ₀, φi, Ny) - 1

    j₁⁻ = floor(Int, j₁)
    j₁⁺ = j₁⁻ + 1
    
    @inbounds begin
        λ₁⁻ = λ[1, j₁⁻]
        λ₁⁺ = λ[1, j₁⁺]
        φ₁⁻ = φ[1, j₁⁻]
        φ₁⁺ = φ[1, j₁⁺]

        for i in 2:Nx 
            # Find j-indices corresponding to φ₀
            φi = view(φ, i, :)
            j₂ = fractional_index(φ₀, φi, Ny) 
            
            j₂⁻ = floor(Int, j₂)
            j₂⁺ = j₂⁻ + 1

            if j₂⁻ > grid.Ny && j₂⁺ > grid.Ny
            
            else
                φ₂⁺ = φ[i, j₂⁺]
                φ₂⁻ = φ[i, j₂⁻]

                # Define the additional 2 points
                λ₂⁺ = λ[i, j₂⁺]
                λ₂⁻ = λ[i, j₂⁻]
                
                # We start by defining the lines y = m x + q that define the rectangle
                # We have to remember that λ could jump between 0 and 360, 
                # so we need to correct for it?

                # Check whether our point is contained in the rectangle
                # described by p₁₁, p₁₂, p₂₁, p₂₂. We orient the points clockwise
                # p₁₁ -> p₁₂ -> p₂₂ -> p₂₁ 

                # Vertical line between p₁₁ and p₁₂
                mᵥ₁ = ifelse(φ₁⁺ == φ₁⁻, zero(φ₁⁻), (λ₁⁺ - λ₁⁻) / (φ₁⁺ - φ₁⁻))
                qᵥ₁ = λ₁⁻ - mᵥ₁ * φ₁⁻

                # Vertical line between p₂₁ and p₂₂
                mᵥ₂ = ifelse(φ₂⁺ == φ₂⁻, zero(φ₂⁻), (λ₂⁺ - λ₂⁻) / (φ₂⁺ - φ₂⁻))
                qᵥ₂ = λ₂⁻ - mᵥ₂ * φ₂⁻

                # vertical bounding lines for λ₀
                λᵥ₁ = mᵥ₁ * φ₀ + qᵥ₁
                λᵥ₂ = mᵥ₂ * φ₀ + qᵥ₂

                # If λᵥ₁ > λᵥ₂ it means we are crossing the 2π line, so
                # we need to correct λᵥ₁ by a factor 2π. 
                # NOTE: this assumes that `OSSG` has coordinates that are expressed
                # in degrees and not in radians.
                λᵥ₁ = ifelse(λᵥ₁ > λᵥ₂, λᵥ₁ - 360, λᵥ₁)

                # @show λᵥ₁, λᵥ₂, λ₀
                # Check that λ₀ lies inbetween the two vertical lines
                horizontal_check = λᵥ₁ ≤ λ₀ ≤ λᵥ₂

                if horizontal_check
                    @show λᵥ₁, λ₀, λᵥ₂
                end
                
                iⁿ = 1 / (λᵥ₂ - λᵥ₁) * (λ₀ - λᵥ₁) + i - 1
                iⁿ = ifelse(λᵥ₂ == λᵥ₁, i, iⁿ)

                # Horizontal line between p₁₁ and p₂₁
                mₕ₁ = (φ₂⁻ - φ₁⁻) / (λ₂⁻ - λ₁⁻) 
                qₕ₁ = φ₁⁻ - mₕ₁ * λ₁⁻

                # Horizontal line between p₁₂ and p₂₂
                mₕ₂ = (φ₂⁺ - φ₁⁺) / (λ₂⁺ - λ₁⁺) 
                qₕ₂ = φ₁⁺ - mₕ₂ * λ₁⁺

                # Horizontal bounding lines for φ₀
                φₕ₁ = mₕ₁ * λ₀ + qₕ₁
                φₕ₂ = mₕ₂ * λ₀ + qₕ₂

                # @show φₕ₁, φₕ₂, φ₀
                # Check that φ₀ lies inbetween the horizontal lines
                vertical_check = φₕ₁ ≤ φ₀ ≤ φₕ₂

                if vertical_check
                    @show φₕ₁, φₕ₂, mₕ₁, mₕ₂
                end

                # Interpolating to find fractional index
                jⁿ⁻ = (j₂⁻ - j₁⁻) / (φ₂⁻ - φ₁⁻) * (φₕ₁ - φ₁⁻) + j₁⁻
                jⁿ⁺ = (j₂⁺ - j₁⁺) / (φ₂⁺ - φ₁⁺) * (φₕ₂ - φ₁⁺) + j₁⁺
                
                jⁿ⁻ = ifelse(φ₂⁻ == φ₁⁻, j₁⁻, jⁿ⁻)
                jⁿ⁺ = ifelse(φ₂⁺ == φ₁⁺, j₁⁺, jⁿ⁺)

                # Final jⁿ index
                jⁿ  = (jⁿ⁺ - jⁿ⁻) / (φₕ₂ - φₕ₁) * (φ₀ - φₕ₁) + jⁿ⁻
                jⁿ  = ifelse(φₕ₂ == φₕ₁, jⁿ⁻, jⁿ)

                update_indices = horizontal_check & vertical_check

                ii = ifelse(update_indices, iⁿ, one(iⁿ))
                jj = ifelse(update_indices, jⁿ, one(jⁿ))
                
                # if update_indices
                #     @show i, jⁿ, j₂⁻, j₂⁺, j₂, φ₂⁺, φ₂⁻, λ₂⁻, λ₂⁺, λ₁⁻, λ₁⁺, λᵥ₁, λᵥ₂
                # end

                # Update indices and
                # previous coordinates
                j₁⁻ = j₂⁻
                j₁⁺ = j₂⁺

                λ₁⁻ = λ₂⁻
                λ₁⁺ = λ₂⁺

                φ₁⁻ = φ₂⁻
                φ₁⁺ = φ₂⁺
            end
        end
    end
   
    return ii, jj
end