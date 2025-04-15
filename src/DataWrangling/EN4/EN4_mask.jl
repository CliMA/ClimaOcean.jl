using Oceananigans: location
using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.Grids: znode

import ClimaOcean: stateindex

"""
    EN4_mask(architecture = CPU(); minimum_value = Float32(-1e5))

A boolean field where `true` represents a missing value in the EN4 dataset.
"""
function EN4_mask(metadata, architecture = CPU();
                   data_field = EN4_field(metadata; architecture, inpainting=nothing),
                   minimum_value = Float32(-1e5),
                   maximum_value = Float32(1e5))

    mask  = Field{location(data_field)...}(data_field.grid, Bool)

    # Set the mask with zeros where field is defined
    launch!(architecture, data_field.grid, :xyz, _set_EN4_mask!, mask, data_field, minimum_value, maximum_value)

    return mask
end

# Default
EN4_mask(arch::AbstractArchitecture=CPU()) = EN4_mask(Metadata(:temperature, dataset=EN4Monthly()), arch)

@kernel function _set_EN4_mask!(mask, Tᵢ, args...)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] == 1e10)
end

"""
    EN4_immersed_grid(metadata, architecture = CPU())

Compute the `ImmersedBoundaryGrid` for `metadata` with a bottom height field that is defined
by the first non-missing value from the bottom up.
"""
function EN4_immersed_grid(metadata, architecture = CPU())

    mask = EN4_mask(metadata, architecture)
    grid = mask.grid
    bottom = Field{Center, Center, Nothing}(grid)

    # Set the mask with zeros where field is defined
    launch!(architecture, grid, :xy, _set_height_from_mask!, bottom, grid, mask)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))
end

# Default
EN4_immersed_grid(arch::AbstractArchitecture=CPU()) = EN4_immersed_grid(Metadata(:temperature, dataset=EN4Monthly()), arch)
