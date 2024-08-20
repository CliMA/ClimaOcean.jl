using OrthogonalSphericalShellGrids
using OrthogonalSphericalShellGrids: TRG
using Oceananigans
using NearestNeighbors
using Oceananigans.Grids: OSSG, λnodes, φnodes
using Oceananigans.Fields: fractional_index, fractional_z_index

import Oceananigans.Fields: interpolate, interpolate!, fractional_indices

ZFlatOSSG = OrthogonalSphericalShellGrid{<:Any, <:Any, <:Any, <:Flat}

@inline function interpolate(at_node, from_field, from_loc, from_grid::OSSG)
    λ₀, φ₀, z₀ = at_node
    i₀, j₀, d₀₀, d₀₁, d₁₀, d₀₂, d₂₀, d₁₁, d₂₂, d₁₂, d₂₁ = horizontal_distances(λ₀, φ₀, from_loc, grid)

    i₁ = i₀ - 1
    i₂ = i₀ + 1

    j₁ = j₀ - 1
    j₂ = j₀ + 1

    @inbounds begin
        f₀₀ = from_field[i₀, j₀, 1]
        f₀₁ = from_field[i₀, j₁, 1]
        f₁₀ = from_field[i₁, j₀, 1]
        f₀₂ = from_field[i₀, j₂, 1]
        f₂₀ = from_field[i₂, j₀, 1]
        f₁₁ = from_field[i₁, j₁, 1]
        f₂₂ = from_field[i₂, j₂, 1]
        f₁₂ = from_field[i₁, j₂, 1]
        f₂₁ = from_field[i₂, j₁, 1]
    end

    w₀₀ = 1 / d₀₀
    w₀₁ = 1 / d₀₁
    w₁₀ = 1 / d₁₀
    w₀₂ = 1 / d₀₂
    w₂₀ = 1 / d₂₀
    w₁₁ = 1 / d₁₁
    w₂₂ = 1 / d₂₂
    w₁₂ = 1 / d₁₂
    w₂₁ = 1 / d₂₁

    f = f₀₀ * w₀₀ + f₀₁ * w₀₁ + f₁₀ * w₁₀ + f₀₂ * w₀₂ + f₂₀ * w₂₀ + f₁₁ * w₁₁ + f₂₂ * w₂₂ + f₁₂ * w₁₂ + f₂₁ * w₂₁
    
    return f / (w₀₀ + w₀₁ + w₁₀ + w₀₂ + w₂₀ + w₁₁ + w₂₂ + w₁₂ + w₂₁)
end

@inline function distance(x₁, x₂, y₁, y₂) 
    dx = x₁ - x₂
    dy = y₁ - y₂
    return dx * dx + dy * dy
end

@inline function check_and_update(dist, i₀, j₀, i, j, λ₀, φ₀, λ, φ)               
    d = distance(λ₀, λ, φ₀, φ) 
    i₀ = ifelse(d < dist, i, i₀)          
    j₀ = ifelse(d < dist, j, j₀)          
    dist = min(d, dist)

    return dist, i₀, j₀
end

# # We assume that in an OSSG, the latitude lines for a given i - index are sorted
# # i.e. φ is monotone in j. This is not the case for λ that might jump between 0 and 360.
@inline function horizontal_distances(λ₀, φ₀, loc, grid)
    # This is a "naive" algorithm, so it is going to be quite slow!
    # Optimizations are welcome!
    λ = λnodes(grid, loc...)
    φ = φnodes(grid, loc...)

    Nx, Ny, _ = size(grid)

    # We search for an initial valid option
    dist = Inf
    i₀ = 1
    j₀ = 1

    @inbounds begin
        for i = 1:Nx
            φi = view(φ, i, :)
            jⁿ = fractional_index(φ₀, φi, Ny) - 1

            j⁻ = floor(Int, jⁿ)
            j⁺ = j⁻ + 1

            if j⁻ <= grid.Ny
                λ⁻ = λ[i, j⁻]
                φ⁻ = φ[i, j⁻]
                dist, i₀, j₀ = check_and_update(dist, i₀, j₀, i, j⁻, λ₀, φ₀, λ⁻, φ⁻)               
            end

            if j⁺ <= grid.Ny
                λ⁺ = λ[i, j⁺]
                φ⁺ = φ[i, j⁺]
                dist, i₀, j₀ = check_and_update(dist, i₀, j₀, i, j⁺, λ₀, φ₀, λ⁺, φ⁺)
            end
        end 
    end
    
    di₀j₀ = dist
    
    # Now find the closest neighbors given i₀ and j₀
    i₁ = i₀ - 1
    j₁ = j₀ - 1
    i₂ = i₀ + 1
    j₂ = j₀ + 1

    di₀j₁ = @inbounds distance(λ₀, φ₀, λ[i₀, j₁], φ[i₀, j₁])
    di₁j₀ = @inbounds distance(λ₀, φ₀, λ[i₁, j₀], φ[i₁, j₀])
    di₀j₂ = @inbounds distance(λ₀, φ₀, λ[i₀, j₂], φ[i₀, j₂])
    di₂j₀ = @inbounds distance(λ₀, φ₀, λ[i₂, j₀], φ[i₂, j₀])
    
    di₁j₁ = @inbounds distance(λ₀, φ₀, λ[i₁, j₁], φ[i₁, j₁])
    di₂j₂ = @inbounds distance(λ₀, φ₀, λ[i₂, j₂], φ[i₂, j₂])
    di₁j₂ = @inbounds distance(λ₀, φ₀, λ[i₁, j₂], φ[i₁, j₂])
    di₂j₁ = @inbounds distance(λ₀, φ₀, λ[i₂, j₁], φ[i₂, j₁])

    return i₀, j₀, di₀j₀, di₀j₁, di₁j₀, di₀j₂, di₂j₀, di₁j₁, di₂j₂, di₁j₂, di₂j₁
end

# # We assume that in an OSSG, the latitude lines for a given i - index are sorted
# # i.e. φ is monotone in j. This is not the case for λ that might jump between 0 and 360.
# @inline function fractional_horizontal_indices(λ₀, φ₀, loc, grid)
#     # This is a "naive" algorithm, so it is going to be quite slow!
#     # Optimizations are welcome!
#     λ = λnodes(grid, loc...)
#     φ = φnodes(grid, loc...)

#     Nx, Ny, _ = size(grid)

#     # Initial indices
#     ii = one(λ₀)
#     jj = one(φ₀)

#     # We search for an initial valid option
#     i₀ = 1 
#     φi = view(φ, i₀, :)
#     j₁ = fractional_index(φ₀, φi, Ny) - 1

#     j₁⁻ = floor(Int, j₁)
#     j₁⁺ = j₁⁻ + 1

#     while j₁⁻ > grid.Ny || j₁⁺ > grid.Ny
#         i₀ += 1 
#         φi = view(φ, i₀, :)
#         j₁ = fractional_index(φ₀, φi, Ny) - 1

#         j₁⁻ = floor(Int, j₁)
#         j₁⁺ = j₁⁻ + 1
#     end

#     @inbounds begin
#         λ₁⁻ = λ[i₀, j₁⁻]
#         λ₁⁺ = λ[i₀, j₁⁺]
#         φ₁⁻ = φ[i₀, j₁⁻]
#         φ₁⁺ = φ[i₀, j₁⁺]

#         # Starting the new loop
#         i = i₀

#         while i <= Nx
#             # Find j-indices corresponding to φ₀
#             φi = view(φ, i, :)
#             j₂ = fractional_index(φ₀, φi, Ny) 
            
#             j₂⁻ = floor(Int, j₂)
#             j₂⁺ = j₂⁻ + 1

#             # Search for the next valid option!
#             if j₂⁻ > grid.Ny && j₂⁺ > grid.Ny
#                 while j₂⁻ > grid.Ny || j₂⁺ > grid.Ny
#                     i += 1
#                     φi = view(φ, i, :)
#                     j₂ = fractional_index(φ₀, φi, Ny) 
                    
#                     j₂⁻ = floor(Int, j₂)
#                     j₂⁺ = j₂⁻ + 1
#                 end        
#             else
#                 φ₂⁺ = φ[i, j₂⁺]
#                 φ₂⁻ = φ[i, j₂⁻]

#                 # Define the additional 2 points
#                 λ₂⁺ = λ[i, j₂⁺]
#                 λ₂⁻ = λ[i, j₂⁻]
                
#                 # We now have our four points, they should be arranged in a rectangle

#                 # We start by defining the lines y = m x + q that define the rectangle
#                 # We have to remember that λ could jump between 0 and 360, 
#                 # so we need to correct for it?

#                 # Check whether our point is contained in the rectangle
#                 # described by p₁₁, p₁₂, p₂₁, p₂₂. We orient the points clockwise
#                 # p₁₁ -> p₁₂ -> p₂₂ -> p₂₁ 
                                
#                 # Fix λ to make sure we account for periodic boundaries
#                 λ₁⁻ = ifelse(λ₁⁻ > λ₂⁻, λ₁⁻ - 360, λ₁⁻)
#                 λ₁⁺ = ifelse(λ₁⁺ > λ₂⁺, λ₁⁺ - 360, λ₁⁺)

#                 # Vertical line between p₁₁ and p₁₂
#                 mᵥ₁ = ifelse(φ₁⁺ == φ₁⁻, zero(φ₁⁻), (λ₁⁺ - λ₁⁻) / (φ₁⁺ - φ₁⁻))
#                 qᵥ₁ = λ₁⁻ - mᵥ₁ * φ₁⁻

#                 # Vertical line between p₂₁ and p₂₂
#                 mᵥ₂ = ifelse(φ₂⁺ == φ₂⁻, zero(φ₂⁻), (λ₂⁺ - λ₂⁻) / (φ₂⁺ - φ₂⁻))
#                 qᵥ₂ = λ₂⁻ - mᵥ₂ * φ₂⁻

#                 # vertical bounding lines for λ₀
#                 λᵥ₁ = mᵥ₁ * φ₀ + qᵥ₁
#                 λᵥ₂ = mᵥ₂ * φ₀ + qᵥ₂

#                 # If λᵥ₁ > λᵥ₂ it means we are crossing the 2π line, so
#                 # we need to correct λᵥ₁ by a factor 2π. 
#                 # NOTE: this assumes that `OSSG` has coordinates that are expressed
#                 # in degrees and not in radians.
#                 # left_of_date_change = λ₀ < 10 && λᵥ₁ > 300
#                 # λᵥ₁ = ifelse(left_of_date_change, λᵥ₁ - 360, λᵥ₁)
#                 # right_of_date_change = λ₀ > 350 && λᵥ₂ < 10
#                 # λᵥ₂ = ifelse(right_of_date_change, λᵥ₂ + λᵥ₁, λᵥ₂)

#                 # @show λᵥ₁, λᵥ₂, λ₀
#                 # Check that λ₀ lies inbetween the two vertical lines
#                 update_indices = λᵥ₁ ≤ λ₀ ≤ λᵥ₂

#                 # @show i
#                 # if update_indices
#                 #     @show λᵥ₁, λᵥ₂
#                 #     @show λ₂⁺, φ₂⁺
#                 #     @show λ₂⁻, φ₂⁻
#                 #     @show λ₁⁻, φ₁⁻
#                 #     @show λ₁⁺, φ₁⁺
#                 # end
                
#                 iⁿ = 1 / (λᵥ₂ - λᵥ₁) * (λ₀ - λᵥ₁) + i - 1
#                 iⁿ = ifelse(λᵥ₂ == λᵥ₁, i, iⁿ)

#                 # We do not need the vertical check because we 
#                 # know that φ₀ is contained in the indices we have 
#                 # calculated. This is because φ is a monotonic function of j
#                 # and we could use the simple 1D index search on the i-column

#                 # Interpolating to find fractional index
#                 jⁿ⁻ = (j₁⁺ - j₁⁻) / (φ₁⁺ - φ₁⁻) * (φ₀ - φ₁⁻) + j₁⁻
#                 jⁿ⁺ = (j₂⁺ - j₂⁻) / (φ₂⁺ - φ₂⁻) * (φ₀ - φ₂⁻) + j₂⁻
                
#                 jⁿ⁻ = ifelse(φ₁⁺ == φ₁⁻, j₁⁻, jⁿ⁻)
#                 jⁿ⁺ = ifelse(φ₂⁺ == φ₂⁻, j₁⁺, jⁿ⁺)

#                 # Final jⁿ index, weigthed on the i-direction
#                 iⁿ⁻ = floor(Int, iⁿ)
#                 jⁿ  = (jⁿ⁺ - jⁿ⁻) * (iⁿ - iⁿ⁻) + jⁿ⁻

#                 ii = ifelse(update_indices, iⁿ, ii)
#                 jj = ifelse(update_indices, jⁿ, jj)

#                 # Update counter
#                 i += 1

#                 # Update indices and
#                 # previous coordinates
#                 j₁⁻ = j₂⁻
#                 j₁⁺ = j₂⁺

#                 λ₁⁻ = λ₂⁻
#                 λ₁⁺ = λ₂⁺

#                 φ₁⁻ = φ₂⁻
#                 φ₁⁺ = φ₂⁺
#             end
#         end
#     end
   
#     return ii, jj
# end