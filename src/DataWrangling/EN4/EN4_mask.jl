using Oceananigans.Architectures: AbstractArchitecture
using ClimaOcean.DataWrangling: dataset_mask, dataset_immersed_grid

import ClimaOcean.DataWrangling: default_set_dataset_mask

# Defaults
EN4_mask(arch::AbstractArchitecture=CPU()) = dataset_mask(Metadatum(:temperature, dataset=EN4Monthly()), arch)
EN4_immersed_grid(arch::AbstractArchitecture=CPU()) = dataset_immersed_grid(Metadatum(:temperature, dataset=EN4Monthly()), arch)

default_set_dataset_mask(metadata::Metadata{<:EN4Monthly}) = ClimaOcean.DataWrangling.EN4._set_EN4_mask!

@kernel function _set_EN4_mask!(mask, Tᵢ, args...)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] == 1e10)
end
