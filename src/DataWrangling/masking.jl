using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Architectures: architecture, device, GPU
using Oceananigans.Fields: Field, OneField, instantiated_location, interior, CenterField
using Oceananigans.Grids: peripheral_node, znode
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

import ClimaOcean: stateindex

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

"""
    NearestNeighborInpainting{M}

A structure representing the nearest neighbor inpainting algorithm, where a missing value is
substituted with the average of the surrounding valid values. This process is repeated a maximum
of `maxiter` times or until the field is completely inpainted.
"""
struct NearestNeighborInpainting{M}
    maxiter :: M
end

# Maybe we can remove this propagate field in lieu of a diffusion,
# Still we'll need to do this a couple of steps on the original grid
@kernel function _propagate_field!(field, ::NearestNeighborInpainting, tmp_field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        nw = field[i - 1, j, k]
        ns = field[i, j - 1, k]
        ne = field[i + 1, j, k]
        nn = field[i, j + 1, k]
        nb = (nw, ne, nn, ns)
    end

    counter = 0
    cumsum  = zero(eltype(field))

    for n in nb
        counter += ifelse(isnan(n), 0, 1)
        cumsum  += ifelse(isnan(n), 0, n)
    end

    @inbounds tmp_field[i, j, k] = ifelse(cumsum == 0, NaN, cumsum / counter)
end

@kernel function _substitute_values!(field, tmp_field)
    i, j, k = @index(Global, NTuple)
    @inbounds needs_inpainting = isnan(field[i, j, k])
    @inbounds field[i, j, k] = ifelse(needs_inpainting, tmp_field[i, j, k], field[i, j, k])
end

@kernel function _nan_mask!(field, mask)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = ifelse(mask[i, j, k], NaN, field[i, j, k])
end

propagate_horizontally!(field, ::Nothing, tmp_field=deepcopy(field); kw...) = field

function propagating(field, mask, iter, inpainting::NearestNeighborInpainting)
    mask_sum = sum(field; condition=interior(mask))
    return isnan(mask_sum) && iter < inpainting.maxiter
end

@kernel function _fill_nans!(field)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = ifelse(isnan(field[i, j, k]), 0, field[i, j, k])
end

"""
    propagate_horizontally!(inpainting, field, mask [, tmp_field=deepcopy(field)])

Horizontally propagate the values of `field` into the `mask`.
In other words, cells where `mask[i, j, k] == false` are preserved,
and cells where `mask[i, j, k] == true` are painted over.

The first argument `inpainting` is the inpainting algorithm to use in the `_propagate_field!` step.
"""
function propagate_horizontally!(inpainting::NearestNeighborInpainting, field, mask, tmp_field=deepcopy(field))
    iter  = 0
    grid  = field.grid
    arch  = architecture(grid)

    launch!(arch, grid, size(field), _nan_mask!, field, mask)
    fill_halo_regions!(field)

    # Need temporary field to avoid a race condition
    parent(tmp_field) .= parent(field)

    while propagating(field, mask, iter, inpainting)
        launch!(arch, grid, size(field), _propagate_field!,   field, inpainting, tmp_field)
        launch!(arch, grid, size(field), _substitute_values!, field, tmp_field)
        iter += 1
        @debug "Propagate pass $iter with sum $(sum(parent(field)))"
    end

    launch!(arch, grid, size(field), _fill_nans!, field)

    fill_halo_regions!(field)

    return field
end

continue_downwards!(field, ::Nothing) = field

"""
    continue_downwards!(field, mask)

Continue downwards a field with missing values within `mask`.
Cells where `mask[i, k, k] == false` will be preserved.
"""
function continue_downwards!(field, mask)
    arch = architecture(field)
    grid = field.grid
    launch!(arch, grid, :xy, _continue_downwards!, field, grid, mask)
    return field
end

@kernel function _continue_downwards!(field, grid, mask)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    for k = Nz-1 : -1 : 1
        @inbounds field[i, j, k] = ifelse(mask[i, j, k], field[i, j, k+1], field[i, j, k])
    end
end

"""
    inpaint_mask!(field, mask; inpainting=NearestNeighborInpainting(Inf))

Inpaint `field` within `mask`, using values outside `mask`.
In other words, regions where `mask[i, j, k] == 1` is inpainted
and regions where `mask[i, j, k] == 0` are preserved.

Arguments
=========

- `field`: `Field` to be inpainted.
- `mask`: Boolean-valued `Field`, values where
          `mask[i, j, k] == true` are inpainted.
- `inpainting`: The inpainting algorithm to use. The only option is
                `NearestNeighborInpainting(maxiter)`, where an average
                of the valid surrounding values is used `maxiter` times.
                Default: `NearestNeighborInpainting(Inf)`.
"""
function inpaint_mask!(field, mask; inpainting=NearestNeighborInpainting(Inf))

    if inpainting isa Int
        inpainting = NearestNeighborInpainting(inpainting)
    end

    # TODO: make this a proper condition
    Nx, Ny, Nz = size(field)
    if Nz > 1
        continue_downwards!(field, mask)
    end

    propagate_horizontally!(inpainting, field, mask)

    return field
end

default_mask_value(metadata) = NaN

"""
    compute_mask(metadata::Metadatum, dataset_field,
                 mask_value = default_mask_value(metadata),
                 minimum_value = -1f5,
                 maximum_value = 1f5)

A boolean field where `true` represents a missing value in the dataset_field.
"""
function compute_mask(metadata::Metadatum, dataset_field,
                      mask_value = default_mask_value(metadata),
                      minimum_value = -1f5,
                      maximum_value = 1f5)

    grid = dataset_field.grid
    arch = Oceananigans.Architectures.architecture(grid)
    LX, LY, LZ = location(dataset_field)
    mask = Field{LX, LY, LZ}(grid, Bool)

    # Set the mask with zeros where field is defined
    launch!(arch, grid, :xyz, _compute_mask!,
            mask, dataset_field, minimum_value, maximum_value, mask_value)

    return mask
end

@kernel function _compute_mask!(mask, field, min_value, max_value, mask_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = is_masked(field[i, j, k], min_value, max_value, mask_value)
end

@inline is_masked(a, min_value, max_value, mask_value) =
    isnan(a) | (a <= min_value) | (a >= max_value) | (a == mask_value)
