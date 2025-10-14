using KernelAbstractions: @kernel, @index

"""
    NearestNeighborInpainting{M}

A structure representing the nearest neighbor inpainting algorithm, where a missing value is
substituted with the average of the surrounding valid values. This process is repeated a maximum
of `maxiter` times or until the field is completely inpainted.
"""
struct NearestNeighborInpainting{M}
    maxiter :: M
end

"""
    ValueInpainting{V}

A structure representing a simple value inpainting algorithm, where all missing values are
substituted with a specified constant value.

# Fields
- `value :: V`: The constant value to use for inpainting missing data.
- `maxiter`: For cache compatibility (not used in the algorithm itself). Defaults to 0.

# Example
```julia
inpainting = ValueInpainting(0.0)  # Fill missing values with 0
inpaint_mask!(field, mask; inpainting)
```
"""
struct ValueInpainting{V, M}
    value :: V
    maxiter :: M
end

# Convenience constructor that sets maxiter to a default value
ValueInpainting(value) = ValueInpainting(value, 0)

propagate_horizontally!(field, ::Nothing, substituting_field=deepcopy(field); kw...) = field

function propagating(field, mask, iter, inpainting::NearestNeighborInpainting)
    nans = sum(isnan, field; condition=interior(mask))
    return nans > 0 && iter < inpainting.maxiter
end

# ValueInpainting doesn't need iteration
propagating(field, mask, iter, inpainting::ValueInpainting) = false

"""
    propagate_horizontally!(inpainting, field, mask [, substituting_field=deepcopy(field)])

Horizontally propagate the values of `field` into the `mask`.
In other words, cells where `mask[i, j, k] == false` are preserved,
and cells where `mask[i, j, k] == true` are painted over.

The first argument `inpainting` is the inpainting algorithm to use in the `_propagate_field!` step.
"""
function propagate_horizontally!(inpainting::NearestNeighborInpainting, field, mask,
                                 substituting_field=deepcopy(field))
    iter  = 0
    grid  = field.grid
    arch  = architecture(grid)

    launch!(arch, grid, size(field), _nan_mask!, field, mask)
    fill_halo_regions!(field)

    # Need temporary field to avoid a race condition
    parent(substituting_field) .= parent(field)

    while propagating(field, mask, iter, inpainting)
        launch!(arch, grid, size(field), _propagate_field!, substituting_field, inpainting, field)
        launch!(arch, grid, size(field), _substitute_values!, field, substituting_field)

        @debug begin
            nans = sum(isnan, field; condition=interior(mask))
            "Propagate pass: $iter, remaining NaNs: $nans"
        end

        iter += 1
    end

    launch!(arch, grid, size(field), _fill_nans!, field)
    fill_halo_regions!(field)

    return field
end

function propagate_horizontally!(inpainting::ValueInpainting, field, mask,
                                 substituting_field=deepcopy(field))
    grid = field.grid
    arch = architecture(grid)

    launch!(arch, grid, size(field), _fill_with_number!, field, mask, inpainting.value)
    fill_halo_regions!(field)

    return field
end

# Maybe we can remove this propagate field in lieu of a diffusion,
# Still we'll need to do this a couple of steps on the original grid
@kernel function _propagate_field!(substituting_field, ::NearestNeighborInpainting, field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        nw = field[i - 1, j, k]
        ns = field[i, j - 1, k]
        ne = field[i + 1, j, k]
        nn = field[i, j + 1, k]
    end

    neighbors = (nw, ne, nn, ns)
    FT = eltype(field)
    donors = 0
    value = zero(FT)

    for n in neighbors
        donors += !isnan(n)
        value += !isnan(n) * n
    end

    FT_NaN = convert(FT, NaN)
    @inbounds substituting_field[i, j, k] = ifelse(value == 0, FT_NaN, value / donors)
end

@kernel function _substitute_values!(field, substituting_field)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        needs_inpainting = isnan(field[i, j, k])
        field[i, j, k] = ifelse(needs_inpainting, substituting_field[i, j, k], field[i, j, k])
    end
end

@kernel function _nan_mask!(field, mask)
    i, j, k = @index(Global, NTuple)
    FT_NaN = convert(eltype(field), NaN)
    @inbounds field[i, j, k] = ifelse(mask[i, j, k], FT_NaN, field[i, j, k])
end

@kernel function _fill_nans!(field)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] *= !isnan(field[i, j, k])
end

@kernel function _fill_with_number!(field, mask, value)
    i, j, k = @index(Global, NTuple)
    FT_value = convert(eltype(field), value)
    @inbounds field[i, j, k] = ifelse(mask[i, j, k], FT_value, field[i, j, k])
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

- `inpainting`: The inpainting algorithm to use. Options include:
                * `NearestNeighborInpainting(maxiter)`: Uses an average
                  of the valid surrounding values, repeated `maxiter` times.
                  Default: `NearestNeighborInpainting(Inf)`.
                * `ValueInpainting(value)`: Fills all missing values with
                  the specified constant `value`.

# Examples
```julia
# Use nearest neighbor inpainting with default settings
inpaint_mask!(field, mask)

# Fill missing values with zero
inpaint_mask!(field, mask; inpainting=ValueInpainting(0))

# Fill missing values with a specific temperature
inpaint_mask!(field, mask; inpainting=ValueInpainting(273.15))
```
"""
function inpaint_mask!(field, mask; inpainting=NearestNeighborInpainting(Inf))

    if inpainting isa Int
        inpainting = NearestNeighborInpainting(inpainting)
    end

    if size(field, 3) > 1
        continue_downwards!(field, mask)
    end

    propagate_horizontally!(inpainting, field, mask)

    return field
end

#####
##### Vertical continuation of fields
#####

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
