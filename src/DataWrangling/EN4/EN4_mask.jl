using Oceananigans: location
using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.Grids: znode
using ClimaOcean.DataWrangling: dataset_mask, _set_height_from_mask!

import ClimaOcean: stateindex

import ClimaOcean.DataWrangling: default_set_dataset_mask

# Default
EN4_mask(arch::AbstractArchitecture=CPU()) = dataset_mask(Metadata(:temperature, dataset=EN4Monthly()), arch)

default_set_dataset_mask(metadata::Metadata{<:EN4Monthly}) = ClimaOcean.DataWrangling.EN4._set_EN4_mask!

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

    mask = dataset_mask(metadata, architecture)
    grid = mask.grid
    bottom = Field{Center, Center, Nothing}(grid)

    # Set the mask with zeros where field is defined
    launch!(architecture, grid, :xy, _set_height_from_mask!, bottom, grid, mask)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))
end

# Default
EN4_immersed_grid(arch::AbstractArchitecture=CPU()) = EN4_immersed_grid(Metadata(:temperature, dataset=EN4Monthly()), arch)
