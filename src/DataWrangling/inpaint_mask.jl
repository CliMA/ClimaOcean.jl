using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields: OneField
using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior, CenterField
using Oceananigans.Architectures: architecture, device, GPU

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
    
    launch!(arch, grid, :xyz, _nan_mask!, field, mask)
    fill_halo_regions!(field)

    # Need temporary field to avoid a race condition
    parent(tmp_field) .= parent(field)

    while propagating(field, mask, iter, inpainting)
        launch!(arch, grid, :xyz, _propagate_field!,   field, inpainting, tmp_field)
        launch!(arch, grid, :xyz, _substitute_values!, field, tmp_field)
        iter += 1
        @debug "Propagate pass $iter with sum $(sum(parent(field)))"
    end

    launch!(arch, grid, :xyz, _fill_nans!, field)

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
    inpaint_mask!(field, mask; max_iter = Inf)

Inpaint `field` within `mask`, using values outside `mask`.
In other words, regions where `mask[i, j, k] == 1` is inpainted
and regions where `mask[i, j, k] == 0` are preserved.

Arguments
=========

- `field`: `Field` to be inpainted.
- `mask`: Boolean-valued `Field`, values where
          `mask[i, j, k] == true` are inpainted.
- `inpainting`: The inpainting algorithm to use. For the moment, the only option is `NearestNeighborInpainting(maxiter)`, 
                where an average of the valid surrounding values is used `maxiter` times.
"""
function inpaint_mask!(field, mask; inpainting=NearestNeighborInpainting(Inf))

    if inpainting isa Int
        inpainting = NearestNeighborInpainting(inpainting)
    end

    continue_downwards!(field, mask)
    propagate_horizontally!(inpainting, field, mask)

    return field
end
