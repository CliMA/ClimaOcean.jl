
import Oceananigans.Fields: interior

interior(u::Union{Array, CuArray, SubArray}) = u

conditioned_mean(T::Field, vol, condition = nothing; kwargs...) = 
    conditioned_mean(interior(T), vol; kwargs...)

function conditioned_mean(T, vol, condition = nothing; dims = Colon())
    Tnew = deepcopy(T)

    condition_field!(Tnew, condition)

    mean_active = sum(interior(vol); dims)
    Tmean       = sum(Tnew .* interior(vol); dims)
    
    return Tmean ./ mean_active
end

condition_field!(f::Union{Array, SubArray}, condition=nothing) = nothing

function condition_field!(f, condition = nothing)
    cf = condition_operand(f, condition, 0)
    materialize_condition!(cf);
    return nothing
end


const c = Center()
const f = Face()

@inline z_bottom(i, j, grid) = znode(i, j, 1, grid, c, c, f)
@inline z_bottom(i, j, grid::ImmersedBoundaryGrid) = @inbounds grid.immersed_boundary.bottom_height[i, j]

@inline bottom(i, j, grid) = znode(i, j, 1, grid, c, c, f)

#####
##### MixedLayerDepthField
#####

# b can be temperature (T) or density (ρ)
@kernel function compute_mld!(h, grid, b, Δb, density_criterion)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz

    b_surface = @inbounds (b[i, j, Nz-2] + b[i, j, Nz-1] + b[i, j, Nz]) / 3

    k_start = Nz - 3
    z_ij    = znode(i, j, k_start, grid, c, c, f)

    @unroll for k in k_start : -1 : 1 # scroll from point just below surface

        b⁺ = @inbounds b[i, j, k+1]
        bᵏ = @inbounds b[i, j, k]

        # If temperature decreases downwards, both are > 0
        # If density increases downwards, both are < 0
        Δb⁺ = b_surface - b⁺
        Δbᵏ = b_surface - bᵏ

        zᵏ = znode(i, j, k, grid, c, c, c)
        Δz⁺ = Δzᶜᶜᶠ(i, j, k+1, grid)

        # Assuming temperature decreases downwards and density increases upwards
        if density_criterion # density criterion
            inside_mixed_layer = (Δb⁺ > - Δb) & (Δbᵏ > - Δb)
            just_below_mixed_layer = (Δb⁺ > - Δb) & (Δbᵏ <= - Δb)
            new_z_ij = zᵏ + (- Δb - Δbᵏ) / (Δb⁺ - Δbᵏ) * Δz⁺
        else # temperature criterion
            inside_mixed_layer = (Δb⁺ < Δb) & (Δbᵏ < Δb)
            just_below_mixed_layer = (Δb⁺ < Δb) & (Δbᵏ >= Δb)
            new_z_ij = zᵏ + (Δb - Δbᵏ) / (Δb⁺ - Δbᵏ) * Δz⁺
        end

        # Linearly interpolate to find mixed layer height

        # Replace z_ij if we found a new mixed layer depth
        replace_z = (just_below_mixed_layer | inside_mixed_layer) & !inactive_node(i, j, k, grid, c, c, c)
        z_ij = ifelse(replace_z, new_z_ij, z_ij)
        if just_below_mixed_layer
            break
        end
    end

    # Note "-" since `h` is supposed to be "depth" rather than "height"
    @inbounds h[i, j, 1] = ifelse(inactive_node(i, j, k_start, grid, c, c, c), zero(grid), - z_ij)
end

struct MixedLayerDepthOperand{B, FT}
    buoyancy_operation :: B
    mixed_layer_buoyancy_differential :: FT
    density_criterion :: Bool
end

Base.summary(op::MixedLayerDepthOperand) = "MixedLayerDepthOperand"

function MixedLayerDepth(grid, tracers; Δb = 0.03, density_criterion = false, kw...)
    b_op    = density_criterion ? tracers.ρ : tracers.T
    operand = MixedLayerDepthOperand(b_op, abs(Δb), density_criterion)
    return Field{Center, Center, Nothing}(grid; operand, kw...)
end

const MixedLayerDepthField = Field{Center, Center, Nothing, <:MixedLayerDepthOperand}

function compute!(h::MixedLayerDepthField, time=nothing)
    arch = architecture(h)
    b    = h.operand.buoyancy_operation
    Δb   = h.operand.mixed_layer_buoyancy_differential
    crit = h.operand.density_criterion
    launch!(arch, h.grid, :xy, compute_mld!, h, h.grid, b, Δb, crit)
    return h
end