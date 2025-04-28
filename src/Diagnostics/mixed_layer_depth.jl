using Oceananigans.Grids: static_column_depthᶜᶜᵃ

mutable struct MixedLayerDepthOperand{FT, B}
    buoyancy_perturbation :: B
    difference_criterion :: FT
end

Base.summary(mldo::MixedLayerDepthOperand) = "MixedLayerDepthOperand"

function MixedLayerDepthOperand(bm, grid, tracers; difference_criterion=1e-4)
    buoyancy_perturbation = buoyancy(bm, grid, tracers)
    difference_criterion = convert(eltype(grid), difference_criterion)
    return MixedLayerDepthOperand(buoyancy_perturbation, difference_criterion)
end

const MixedLayerDepthField = Field{<:Any, <:Any, <:Any, <:MixedLayerDepthOperand}

"""
    MixedLayerDepthField(bm, grid, tracers; difference_criterion=1e-4)

"""
function MixedLayerDepthField(bm, grid, tracers; difference_criterion=3e-5)
    operand = MixedLayerDepthOperand(bm, grid, tracers; difference_criterion)
    loc = (Center, Center, Nothing)
    indices = (:, :, :)
    bcs = FieldBoundaryConditions(grid, loc)
    data = new_data(grid, loc, indices)
    recompute_safely = false
    status = FieldStatus()
    return Field(loc, grid, data, bcs, indices, operand, status)
end

function compute!(mld::MixedLayerDepthField, time=nothing)
    compute_mixed_layer_depth!(mld)
    #@apply_regionally compute_mixed_layer_depth!(mld)
    fill_halo_regions!(mld)
    return mld
end

function compute_mixed_layer_depth!(mld)
    grid = mld.grid
    arch = architecture(grid)

    launch!(arch, mld.grid, :xy,
            _compute_mixed_layer_depth!,
            mld,
            grid,
            mld.operand.buoyancy_perturbation,
            mld.operand.difference_criterion)

    return mld
end

const c = Center()
const f = Face()

@kernel function _compute_mixed_layer_depth!(mld, grid, b, Δb★)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    Δb = zero(grid)
    bN = @inbounds b[i, j, Nz]
    mixed = true
    k = Nz - 1
    inactive = inactive_cell(i, j, k, grid)

    while !inactive & mixed & (k > 0)
        Δb = @inbounds bN - b[i, j, k]
        mixed = Δb < Δb★
        k -= 1
        inactive = inactive_cell(i, j, k, grid)
    end

    # Linearly interpolate
    # z★ = zk + Δz/Δb * (Δb★ - Δb)
    zN = znode(i, j, Nz, grid, c, c, c)
    zk = znode(i, j, k, grid, c, c, c)
    Δz = zN - zk
    z★ = zk - Δz/Δb * (Δb★ - Δb)

    # Special case when domain is one grid cell deep
    z★ = ifelse(Δb == 0, zN, z★)

    # Apply various criterion
    h = -z★
    h = max(h, zero(grid))
    H = static_column_depthᶜᶜᵃ(i, j, grid)
    h = min(h, H)

    @inbounds mld[i, j, 1] = h
end

