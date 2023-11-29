using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields: OneField
using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior, CenterField
using Oceananigans.Architectures: architecture, device, GPU

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

# Maybe we can remove this propagate field in lieu of a diffusion, 
# Still we'll need to do this a couple of steps on the original grid
@kernel function _propagate_field!(field, tmp_field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        nw = field[i - 1, j, k]
        ns = field[i, j - 1, k]
        ne = field[i + 1, j, k]
        nn = field[i, j + 1, k]
        nb = (nw, ne, nn, ns)

        counter = 0
        cumsum  = 0.0 

        @unroll for n in nb
            counter += ifelse(isnan(n), 0, 1)
            cumsum  += ifelse(isnan(n), 0, n)
        end

        tmp_field[i, j, k] = ifelse(cumsum == 0, NaN, cumsum / counter)
    end
end

@kernel function _substitute_values!(field, tmp_field)
    i, j, k = @index(Global, NTuple)
    @inbounds substitute = isnan(field[i, j, k])
    @inbounds field[i, j, k] = ifelse(substitute, tmp_field[i, j, k], field[i, j, k])
end

@kernel function _nans_outside_mask!(field, mask)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = ifelse(mask[i, j, k] == 0, NaN, field[i, j, k])
end

propagate_horizontally!(field, ::Nothing; kw...) = nothing

""" 
    propagate_horizontally!(field, mask; max_iter = Inf)

propagate horizontally a field with missing values outside of a `mask`.
Grid cells where `mask == 1` will be preserved
"""
function propagate_horizontally!(field, mask; max_iter = Inf) 
    iter  = 0
    grid  = field.grid
    arch  = architecture(grid)
    
    launch!(arch, grid, :xyz, _nans_outside_mask!, field, mask)
    fill_halo_regions!(field)

    tmp_field = deepcopy(field)

    while isnan(sum(interior(field))) && iter < max_iter
        launch!(arch, grid, :xyz, _propagate_field!,   field, tmp_field)
        launch!(arch, grid, :xyz, _substitute_values!, field, tmp_field)
        iter += 1
        @debug "propagate pass $iter with sum $(sum(parent(field)))"
    end

    return nothing
end

continue_downwards!(field, ::Nothing) = nothing

""" 
    continue_downwards!(field, mask)

continue downwards a field with missing values outside of a `mask`.
Grid cells where `mask == 1` will be preserved
"""
function (field, mask)
    arch = architecture(field)
    grid = field.grid
    launch!(arch, grid, :xy, _continue_downwards!, field, grid, mask)
    return nothing
end

@kernel function _continue_downwards!(field, grid, mask)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz

    @unroll for k = Nz-1 : -1 : 1
        @inbounds fill_from_above = mask[i, j, k] == 0
        @inbounds field[i, j, k] = ifelse(fill_from_above, field[i, j, k+1], field[i, j, k])
    end
end

function fill_missing_values!(tracer; mask = nothing, max_iter = Inf)

    continue_downwards!(tracer, mask)
    propagate_horizontally!(tracer, mask; max_iter)
    
    return tracer
end