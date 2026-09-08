#####
##### Diagnostic restoring of the Denmark Strait overflow product
#####
#
# Not a parameterization. This forces the density delivered to the East Greenland slope to the observed
# value, so that the question "does the overflow deficit set the AMOC?" can be answered without first
# solving "how do we get dense water down a staircase?".
#
# The two bottom boundary layer schemes both left the delivered density at 27.80 against the control's
# 27.80, so neither of them tested the hypothesis — they tested themselves. This does test it, at the
# cost of being unphysical: inside the mask the tracers are simply pinned.
#
# Measured on eORCA1: the mask holds 158 cells, 99 300 km³, spanning 1527-2729 m, currently carrying
# Θ = 3.39, Sᴬ = 35.13, σθ = 27.818. Observed Denmark Strait Overflow Water at this stage of its descent
# is Θ ≈ 2.0, Sᴬ ≈ 35.07, σθ ≈ 27.892 — dense enough to reach ~3800 m rather than the ~1800 m the model
# currently manages.

using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: Center, λnode, φnode, znode
using Oceananigans.ImmersedBoundaries: inactive_node

"""
    overflow_restoring_mask(grid; longitude, latitude, minimum_depth)

Unit mask over the wet cells of the descent region, zero elsewhere. Built on a CPU copy of the grid
because the immersed boundary is not scalar-indexable on a GPU.
"""
function overflow_restoring_mask(grid; longitude, latitude, minimum_depth)
    cpu_grid = on_architecture(CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    values = zeros(eltype(cpu_grid), Nx, Ny, Nz)

    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        inactive_node(i, j, k, cpu_grid, Center(), Center(), Center()) && continue

        λ = λnode(i, j, k, cpu_grid, Center(), Center(), Center())
        λ = ifelse(λ > 180, λ - 360, λ)
        φ = φnode(i, j, k, cpu_grid, Center(), Center(), Center())
        z = znode(i, j, k, cpu_grid, Center(), Center(), Center())

        inside = longitude[1] <= λ <= longitude[2] &&
                 latitude[1]  <= φ <= latitude[2]  &&
                 -z > minimum_depth

        inside && (values[i, j, k] = 1)
    end

    mask = Field{Center, Center, Center}(grid)
    set!(mask, values)
    fill_halo_regions!(mask)

    return mask
end

"""
    overflow_restoring_tendency(i, j, k, grid, clock, fields, parameters)

Discrete-form forcing, `mask * (target - c) / timescale`. Branch-free: outside the mask the rate is
multiplied by zero rather than skipped.
"""
@inline function overflow_restoring_tendency(i, j, k, grid, clock, fields, parameters)
    c = tracer_field(fields, parameters.tracer_name)

    @inbounds begin
        m  = parameters.mask[i, j, k]
        cᵢ = c[i, j, k]
    end

    return m * parameters.rate * (parameters.target - cᵢ)
end

"""
    overflow_restoring_forcing(grid, timescale;
                               longitude, latitude, minimum_depth,
                               target_temperature, target_salinity)

Return `(T, S)` forcings pinning the descent region to observed overflow water, or an empty
`NamedTuple` when `timescale` is `nothing`. Defaults are the observed Denmark Strait Overflow Water
properties and the region measured to carry the model's failed plume.
"""
function overflow_restoring_forcing(grid, timescale;
                                    longitude = (-36.0, -26.0),
                                    latitude = (62.0, 66.5),
                                    minimum_depth = 1500,
                                    target_temperature = 2.0,
                                    target_salinity = 35.065)

    isnothing(timescale) && return NamedTuple()

    mask = overflow_restoring_mask(grid; longitude, latitude, minimum_depth)
    rate = 1 / timescale

    T = Forcing(overflow_restoring_tendency; discrete_form = true,
                parameters = (; mask, rate, target = target_temperature, tracer_name = Val(:T)))

    S = Forcing(overflow_restoring_tendency; discrete_form = true,
                parameters = (; mask, rate, target = target_salinity, tracer_name = Val(:S)))

    return (; T, S)
end
