#####
##### Diagnostic restoring of the Labrador interior salinity toward observation
#####
#
# Not a parameterization. Campaign 27 priced the convection: removing d metres of freshwater from the
# deep Labrador column takes the autumn convective resistance D(1000 m) from `orca`'s 1.170 m² s⁻² to
# `noicedyn`'s 0.750 at d = 1.54-1.65 m, and the model's measured freshwater-content excess over WOA in
# that box (0-200 m, C22-7) is 1.638 m. Those are the same number. This pins the box to the observed
# salinity so that "does the Labrador freshwater bias set the AMOC?" can be answered without first
# solving what causes the bias.
#
# It has to be a restoring rather than a one-off correction because the excess is maintained laterally:
# `orca_davisarch` stopped adding 177 g/kg m of freshwater-equivalent to `labint` west over sixteen
# years and its column kept 2 of it — the basin replaced the rest within a year. A fixed flux would be
# consumed by that compensation; a restoring supplies exactly as much as the compensation removes, and
# the tendency it has to sustain is itself the measurement of how hard the basin fights back.
#
# ⚠ The target is WOA **Annual**, a climatological mean, so inside the box this also damps the
# interannual salinity variability. That is the cost of the attribution test, and it is the reason the
# box is kept to the deep interior (bottom below 2000 m) rather than the whole Labrador: the shelf is
# 62% too fresh for reasons that are not the interior's, and the convecting columns are all interior.

using Oceananigans
using Oceananigans.Architectures: on_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: CenterField, interior
using Oceananigans.Grids: Center, λnode, φnode, znode
using Oceananigans.ImmersedBoundaries: inactive_node
using NumericalEarth.DataWrangling: Metadatum, WOAAnnual

"""
    labrador_restoring_mask(grid; longitude, latitude, minimum_bottom_depth, maximum_depth)

Unit mask over the wet cells of the deep Labrador interior shallower than `maximum_depth`, zero
elsewhere. A column qualifies only if it carries an active cell below `minimum_bottom_depth`, which is
the bathymetric `labint` definition every campaign-27 diagnostic uses, so the restored cells and the
diagnosed cells are the same set. Built on a CPU copy of the grid because the immersed boundary is not
scalar-indexable on a GPU.
"""
function labrador_restoring_mask(grid; longitude, latitude, minimum_bottom_depth, maximum_depth)
    cpu_grid = on_architecture(CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    values = zeros(eltype(cpu_grid), Nx, Ny, Nz)

    for j in 1:Ny, i in 1:Nx
        deep = false
        for k in 1:Nz
            inactive_node(i, j, k, cpu_grid, Center(), Center(), Center()) && continue
            z = znode(i, j, k, cpu_grid, Center(), Center(), Center())
            -z > minimum_bottom_depth && (deep = true; break)
        end
        deep || continue

        λ = λnode(i, j, 1, cpu_grid, Center(), Center(), Center())
        λ = ifelse(λ > 180, λ - 360, λ)
        φ = φnode(i, j, 1, cpu_grid, Center(), Center(), Center())

        (longitude[1] <= λ <= longitude[2] && latitude[1] <= φ <= latitude[2]) || continue

        for k in 1:Nz
            inactive_node(i, j, k, cpu_grid, Center(), Center(), Center()) && continue
            z = znode(i, j, k, cpu_grid, Center(), Center(), Center())
            -z <= maximum_depth && (values[i, j, k] = 1)
        end
    end

    mask = Field{Center, Center, Center}(grid)
    set!(mask, values)
    fill_halo_regions!(mask)

    return mask
end

"""
    labrador_restoring_tendency(i, j, k, grid, clock, fields, parameters)

Discrete-form forcing, `mask * (target - S) / timescale`, with `target` a field rather than a scalar so
the observed vertical structure of the top 200 m is preserved. Pinning the box to a single value would
erase the haline stratification by fiat and guarantee the convection this is meant to test for.
Branch-free: outside the mask the rate is multiplied by zero rather than skipped.
"""
@inline function labrador_restoring_tendency(i, j, k, grid, clock, fields, parameters)
    c = tracer_field(fields, parameters.tracer_name)

    @inbounds begin
        m  = parameters.mask[i, j, k]
        cᵢ = c[i, j, k]
        cᵗ = parameters.target[i, j, k]
    end

    return m * parameters.rate * (cᵗ - cᵢ)
end

"""
    labrador_restoring_forcing(grid, timescale; longitude, latitude, minimum_bottom_depth,
                               maximum_depth, restoring_dir)

Return an `S` forcing pinning the deep Labrador interior above `maximum_depth` to WOA Annual Absolute
Salinity, or an empty `NamedTuple` when `timescale` is `nothing`. Salinity only: campaign 27 measured
the year-1 resistance change as 100% haline (haline-only 0.785 against a full 0.783, thermal 0.998), so
restoring temperature as well would add a term already known to be inert and confound the attribution.

`restoring_dir` must already hold WOA Annual — this runs before the simulation builder creates that
directory, so it reads what the initial condition would read rather than downloading it.
"""
function labrador_restoring_forcing(grid, timescale;
                                    longitude = (-65.0, -40.0),
                                    latitude = (52.0, 66.0),
                                    minimum_bottom_depth = 2000,
                                    maximum_depth = 200,
                                    restoring_dir = "climatology")

    isnothing(timescale) && return NamedTuple()

    mask = labrador_restoring_mask(grid; longitude, latitude, minimum_bottom_depth, maximum_depth)

    # WOA Annual arrives as in-situ temperature and Practical Salinity; the ocean carries Conservative
    # Temperature and Absolute Salinity, and `woa_to_teos10!` converts both in place. The temperature
    # field is needed for the conversion and discarded afterwards.
    T_woa = CenterField(grid)
    target = CenterField(grid)
    set!(T_woa,  Metadatum(:temperature; dir = restoring_dir, dataset = WOAAnnual()))
    set!(target, Metadatum(:salinity;    dir = restoring_dir, dataset = WOAAnnual()))
    woa_to_teos10!(T_woa, target)

    # WOA is not defined in every wet model cell; restoring toward a missing value would inject NaN.
    mask_h   = Array(interior(mask))
    target_h = Array(interior(target))
    mask_h[.!isfinite.(target_h)] .= 0
    copyto!(interior(mask), mask_h)
    fill_halo_regions!(mask)
    fill_halo_regions!(target)

    nwet = count(>(0), mask_h)
    bias = nwet == 0 ? NaN :
           sum(target_h[mask_h .> 0]) / nwet
    @info "Labrador salinity restoring: lon $(longitude), lat $(latitude), bottom below " *
          "$(minimum_bottom_depth) m, above $(maximum_depth) m, timescale $(timescale). " *
          "$(nwet) wet cells restored, WOA mean Sᴬ there = $(round(bias, digits = 4)) g/kg."

    S = Forcing(labrador_restoring_tendency; discrete_form = true,
                parameters = (; mask, target, rate = 1 / timescale, tracer_name = Val(:S)))

    return (; S)
end
