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
    launch!(architecture, data_field.grid, :xyz, _set_mask!, mask, data_field, minimum_value, maximum_value)

    return mask
end

# Default
EN4_mask(arch::AbstractArchitecture=CPU()) = EN4_mask(Metadata(:temperature, dataset=EN4Monthly()), arch)

@kernel function _set_mask!(mask, Tᵢ, args...)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = ismissing(Tᵢ[i, j, k])
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

struct LinearlyTaperedPolarMask{N, S, Z} 
    northern :: N
    southern :: S
    z :: Z
end

"""
    LinearlyTaperedPolarMask(; northern = (70,   75),
                               southern = (-75, -70),
                               z = (-20, 0))

Build a mask that is linearly tapered in latitude between the northern and southern edges.
The mask is constant in depth between the z and equals zero everywhere else.
The mask is limited to lie between (0, 1).
The mask has the following functional form:

```julia
n = 1 / (northern[2] - northern[1]) * (φ - northern[1])
s = 1 / (southern[1] - southern[2]) * (φ - southern[2])

valid_depth = (z[1] < z < z[2])

mask = valid_depth ? clamp(max(n, s), 0, 1) : 0
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
    
    # The mask is active only between `mask.z[1]` and `mask.z[2]`
    valid_depth = (mask.z[1] < z < mask.z[2])

    # we clamp the mask between 0 and 1
    mask_value = clamp(max(n, s), 0, 1)

    return ifelse(valid_depth, mask_value, zero(n))
end

@inline function stateindex(mask::LinearlyTaperedPolarMask, i, j, k, grid, time, loc)
    LX, LY, LZ = loc 
    λ, φ, z = node(i, j, k, grid, LX(), LY(), LZ())
    return mask(φ, z)
end
