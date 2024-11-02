using Oceananigans.Architectures: AbstractArchitecture

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
