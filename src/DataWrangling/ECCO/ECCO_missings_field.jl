using Oceananigans.Architectures: AbstractArchitecture

"""
    ECCO_missings_field([metadata=ECCOMetadata(:temperature),] architecture = CPU();
                        minimum_value = Float32(-1e5),
                        maximum_value = Float32(1e5),
                        filename = metadata_filename(metadata))

A boolean field where `true` represents a missing value in the ECCO dataset.
"""
function ECCO_missings_field(metadata, architecture = CPU(); 
                             minimum_value = Float32(-1e5),
                             maximum_value = Float32(1e5),
                             filename = metadata_filename(metadata))

    field = ecco_field(metadata; architecture, filename)
    mask  = Field{location(field)...}(field.grid, Bool)

    # ECCO4 has zeros in place of the missing values, while
    # ECCO2 expresses missing values with values < -1e5
    if metadata.version isa ECCO4Monthly
        _find_missings! = _find_ecco4_missings!
    else
        _find_missings! = _find_ecco4_missings!
    end

    # Set the mask with zeros where field is defined
    launch!(architecture, field.grid, :xyz,
            _find_missings!, mask, field, minimum_value, maximum_value)

    return mask
end

# Default
ECCO_missings_field(arch::AbstractArchitecture=CPU()) = ECCO_missings_field(ECCOMetadata(:temperature), arch)

@kernel function _find_ecco2_missings!(mask, Tᵢ, minimum_value, maximum_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] < minimum_value) | (Tᵢ[i, j, k] > maximum_value) 
end

@kernel function _find_ecco4_missings!(mask, Tᵢ, minimum_value, maximum_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] == 0) 
end

