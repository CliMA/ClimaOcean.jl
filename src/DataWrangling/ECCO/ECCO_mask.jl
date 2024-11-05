using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.Grids: znode

"""
    ECCO_mask(architecture = CPU(); minimum_value = Float32(-1e5))

A boolean field where `true` represents a missing value in the ECCO dataset.
"""
function ECCO_mask(metadata, architecture = CPU(); 
                   minimum_value = Float32(-1e5),
                   maximum_value = Float32(1e5))

    field = ECCO_field(metadata; architecture)
    mask  = Field{location(field)...}(field.grid, Bool)

    # ECCO4 has zeros in place of the missing values, while
    # ECCO2 expresses missing values with values < -1e5
    if metadata.version isa ECCO4Monthly 
        _set_mask! = _set_ECCO4_mask!
    else
        _set_mask! = _set_ECCO2_mask!
    end

    # Set the mask with zeros where field is defined
    launch!(architecture, field.grid, :xyz, _set_mask!, mask, field, minimum_value, maximum_value)

    return mask
end

# Default
ECCO_mask(arch::AbstractArchitecture=CPU()) = ECCO_mask(ECCOMetadata(:temperature), arch)

@kernel function _set_ECCO2_mask!(mask, Tᵢ, minimum_value, maximum_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] < minimum_value) | (Tᵢ[i, j, k] > maximum_value) 
end

@kernel function _set_ECCO4_mask!(mask, Tᵢ, args...)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] == 0) 
end

"""
    ECCO_immersed_grid(metadata, architecture = CPU())

Compute the `ImmersedBoundaryGrid` for `metadata` with a bottom height field that is defined 
by the first non-missing value from the bottom up.
"""
function ECCO_immersed_grid(metadata, architecture = CPU())

    mask = ECCO_mask(metadata, architecture)
    grid = mask.grid
    bottom = Field{Center, Center, Nothing}(grid)

    # Set the mask with zeros where field is defined
    launch!(architecture, grid, :xy, _set_height_from_mask!, bottom, grid, mask)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))
end

# Default
ECCO_immersed_grid(arch::AbstractArchitecture=CPU()) = ECCO_immersed_grid(ECCOMetadata(:temperature), arch)

@kernel function _set_height_from_mask!(bottom, grid, mask)
    i, j = @index(Global, NTuple)
    
    # Starting from the bottom
    @inbounds bottom[i, j, 1] = znode(i, j, 1, grid, Center(), Center(), Face())

    # Sweep up
    for k in 1:grid.Nz
        z⁺ = znode(i, j, k+1, grid, Center(), Center(), Face())
        @inbounds bottom[i, j, k] = ifelse(mask[i, j, k], z⁺, bottom[i, j, k])
    end
end
