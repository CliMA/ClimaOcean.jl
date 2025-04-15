using Oceananigans: location
using Oceananigans.Architectures: AbstractArchitecture
using Oceananigans.Grids: znode
using ClimaOcean.DataWrangling: dataset_mask, _set_height_from_mask!

import ClimaOcean.DataWrangling: default_set_dataset_mask

# Default
ECCO_mask(arch::AbstractArchitecture=CPU()) = dataset_mask(Metadata(:temperature, dataset=ECCO4Monthly()), arch)

default_set_dataset_mask(metadata::Metadata{<:ECCO4Monthly}) = ClimaOcean.DataWrangling.ECCO._set_ECCO4_mask!
default_set_dataset_mask(metadata::Metadata{<:Union{ECCO2Monthly, ECCO2Daily}}) = ClimaOcean.DataWrangling.ECCO._set_ECCO2_mask!

# ECCO2 expresses missing values with values < -1e5
@kernel function _set_ECCO2_mask!(mask, Tᵢ, minimum_value, maximum_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] < minimum_value) | (Tᵢ[i, j, k] > maximum_value)
end

# ECCO4 has zeros in place of the missing values, while
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

    mask = dataset_mask(metadata, architecture)
    grid = mask.grid
    bottom = Field{Center, Center, Nothing}(grid)

    # Set the mask with zeros where field is defined
    launch!(architecture, grid, :xy, _set_height_from_mask!, bottom, grid, mask)

    return ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))
end

# Default
ECCO_immersed_grid(arch::AbstractArchitecture=CPU()) = ECCO_immersed_grid(Metadata(:temperature, dataset=ECCO4Monthly()), arch)
