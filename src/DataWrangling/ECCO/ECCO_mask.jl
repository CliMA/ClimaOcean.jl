using Oceananigans.Architectures: AbstractArchitecture
using ClimaOcean.DataWrangling: dataset_mask, dataset_immersed_grid

import ClimaOcean.DataWrangling: default_set_dataset_mask

# Defaults
ECCO_mask(arch::AbstractArchitecture=CPU()) = dataset_mask(Metadata(:temperature, dataset=ECCO4Monthly()), arch)
ECCO_immersed_grid(arch::AbstractArchitecture=CPU()) = dataset_immersed_grid(Metadata(:temperature, dataset=ECCO4Monthly()), arch)

default_set_dataset_mask(metadata::Metadata{<:Union{ECCO4Monthly, ECCO4DarwinMonthly}}) = ClimaOcean.DataWrangling.ECCO._set_ECCO4_mask!
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
