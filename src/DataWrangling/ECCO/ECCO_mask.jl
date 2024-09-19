using Oceananigans.Architectures: AbstractArchitecture

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

struct LatitudinallyTaperedPolarMask{N, S, Z} <: Function
    northern_edges :: N
    southern_edges :: S
    z_edges :: F
end

"""
    LatitudinallyTaperedPolarMask(; northern_edges = (70,   75),
                               southern_edges = (-75, -70),
                               z_edges = (-20, 0))

Build a mask that is linearly tapered in latitude between the northern and southern edges.
The mask is constant in depth between the z_edges and is equal to zero everywhere else.
"""
function LatitudinallyTaperedPolarMask(; northern_edges = (70, 75),
                                     southern_edges = (-75, -70),
                                     z_edges = (-20, 0))

    return LatitudinallyTaperedPolarMask(northern_edges, southern_edges, z_edges)
end

@inline function (mask::LinearlyTaperedPolarMask)(λ, φ, z, args...)
    n = 1 / (mask.northern_edges[2] - mask.northern_edges[1]) * (φ - mask.northern_edges[1])
    s = 1 / (mask.southern_edges[1] - mask.southern_edges[2]) * (φ - mask.southern_edges[2])(φ)  
    
    within_depth = (mask.z_edges[1] < z < mask.z_edges[2])

    return ifelse(within_depth, max(n, s, zero(n)), zero(n))
end