using Oceananigans.Architectures: AbstractArchitecture
import ClimaOcean: stateindex

"""
    ECCO_mask(architecture = CPU(); minimum_value = Float32(-1e5))

A boolean field where `true` represents a missing value in the ECCO dataset.
"""
function ECCO_mask(metadata, architecture = CPU(); 
                   minimum_value = Float32(-1e5),
                   maximum_value = Float32(1e5),
                   filename = metadata_filename(metadata))

    field = ECCO_field(metadata; architecture, filename)
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

struct LinearlyTaperedPolarMask{N, S, Z} 
    northern :: N
    southern :: S
    z :: Z
end

"""
    LinearlyTaperedPolarMask(; northern = (70,   75),
                               southern = (-75, -70),
                               z = (-20, 0))

Build a mask that is linearly tapered in latitude inbetween the northern and southern edges.
The mask is constant in depth between the z and is equal to zero everywhere else.
The mask has the following functional form:

```julia
n = 1 / (northern[2] - northern[1]) * (φ - northern[1])
s = 1 / (southern[1] - southern[2]) * (φ - southern[2])

within_depth = (z[1] < z < z[2])

mask = within_depth ? max(n, s, 0) : 0
```
"""
function LinearlyTaperedPolarMask(; northern = (70,   75),
                                    southern = (-75, -70),
                                    z = (-20, 0))

    northern[1] > northern[2]  && throw(ArgumentError("Northern latitude range is invalid, northern[1] > northern[2]."))
    southern[1] > southern[2]  && throw(ArgumentError("Southern latitude range is invalid, southern[1] > southern[2]."))
    z[1] > z[2]                && throw(ArgumentError("Depth range is invalid, z[1] > z[2]."))

    return LinearlyTaperedPolarMask(northern, southern, z)
end

@inline function (mask::LinearlyTaperedPolarMask)(φ, z)
    n = 1 / (mask.northern[2] - mask.northern[1]) * (φ - mask.northern[1])
    s = 1 / (mask.southern[1] - mask.southern[2]) * (φ - mask.southern[2])
    
    within_depth = (mask.z[1] < z < mask.z[2])

    return ifelse(within_depth, max(n, s, zero(n)), zero(n))
end

@inline function stateindex(mask::LinearlyTaperedPolarMask, i, j, k, grid, time, loc)
    LX, LY, LZ = loc 
    λ, φ, z = node(i, j, k, grid, LX(), LY(), LZ())
    return mask(φ, z)
end