#####
##### Diffusive bottom boundary layer (Beckmann & Döscher 1997)
#####
#
# Ported from NEMO `src/OCE/TRA/trabbl.F90` with `nn_bbl_ldf = 1`. Design notes and the measurements
# motivating it are in `docs/plans/2026-08-17-bottom-boundary-layer-design.md`.
#
# The scheme is a two-dimensional Laplacian of the *bottom-cell* tracer values, added to the bottom cell.
# It is not a three-dimensional flux between cells at differing `k`.

using Oceananigans
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.Utils: launch!
using KernelAbstractions: @index, @kernel
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: Center, znode
using Oceananigans.ImmersedBoundaries: inactive_node
using Oceananigans.Operators: Azᶜᶜᶜ, Δxᶜᶠᶜ, Δxᶠᶜᶜ, Δyᶜᶠᶜ, Δyᶠᶜᶜ, Δzᶜᶜᶜ
using SeawaterPolynomials: haline_contraction, thermal_expansion
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Adapt: Adapt, adapt

"""
    BottomBoundaryLayer(grid, equation_of_state; diffusivity)

Dense water sitting upslope of a deeper neighbour is diffused along the bottom, mimicking a gravity
current that a z-coordinate model otherwise mixes away over one or two grid cells.

`bottom_index` is the only precomputed field. Everything else — the slope sign and the buoyancy gradient
that together form the activation criterion — is derived inside the kernel from the same `(i,j)` to
neighbour pairing, so both flip sign together across the tripolar fold and their product is invariant.

`diffusivity` is NEMO's `rn_ahtbbl` in m² s⁻¹. Because this is a symmetric exchange it drives the two
bottom cells to their volume-weighted mean, which caps the density it can deliver downslope regardless of
`diffusivity`; [`AdvectiveBottomBoundaryLayer`](@ref) has no such cap.
"""
struct BottomBoundaryLayer{K, T, FT, E}
    bottom_index      :: K
    transport_x       :: T
    transport_y       :: T
    diffusivity       :: FT
    equation_of_state :: E
end

Adapt.adapt_structure(to, bbl::BottomBoundaryLayer) =
    BottomBoundaryLayer(adapt(to, bbl.bottom_index),
                        adapt(to, bbl.transport_x),
                        adapt(to, bbl.transport_y),
                        adapt(to, bbl.diffusivity),
                        adapt(to, bbl.equation_of_state))

function BottomBoundaryLayer(grid, equation_of_state; diffusivity)

    bottom_index = Field{Center, Center, Nothing}(grid)
    set!(bottom_index, deepest_wet_level(grid))
    fill_halo_regions!(bottom_index)

    transport_x = Field{Center, Center, Nothing}(grid)
    transport_y = Field{Center, Center, Nothing}(grid)

    return BottomBoundaryLayer(bottom_index, transport_x, transport_y,
                               convert(eltype(grid), diffusivity), equation_of_state)
end

"""
    update_bottom_boundary_layer!(sim)

Refresh the two face-transport fields from the current ocean state. Called once per step, so the
activation criterion — which costs four equation-of-state evaluations per face — is paid on the 2D
bottom surface rather than on every cell of the 3D domain for every tracer.
"""
function update_bottom_boundary_layer!(sim, bbl::BottomBoundaryLayer)
    ocean = sim.model.ocean
    fields = (T = ocean.model.tracers.T, S = ocean.model.tracers.S)
    return update_bottom_boundary_layer!(bbl, ocean.model.grid, fields)
end

function update_bottom_boundary_layer!(bbl::BottomBoundaryLayer, grid, fields)
    launch!(architecture(grid), grid, :xy, _compute_bottom_boundary_layer_transport!,
            bbl.transport_x, bbl.transport_y, grid, bbl, fields)

    fill_halo_regions!(bbl.transport_x)
    fill_halo_regions!(bbl.transport_y)

    return nothing
end

@kernel function _compute_bottom_boundary_layer_transport!(Qx, Qy, grid, bbl, fields)
    i, j = @index(Global, NTuple)
    @inbounds begin
        Qx[i, j, 1] = x_face_transport(bbl, fields, grid, i, j)
        Qy[i, j, 1] = y_face_transport(bbl, fields, grid, i, j)
    end
end

# Deepest wet level of every column, zero over land. Built on a CPU copy of the grid because the
# immersed boundary is not scalar-indexable on a GPU.
function deepest_wet_level(grid)
    cpu_grid = on_architecture(CPU(), grid)
    Nx, Ny, Nz = size(cpu_grid)
    levels = zeros(eltype(cpu_grid), Nx, Ny, 1)

    # Nested rather than flattened: a comma-separated `for` shares one `break`, which would stop the sweep
    # after the first wet cell in the whole domain.
    for j in 1:Ny
        for i in 1:Nx
            for k in 1:Nz
                if !inactive_node(i, j, k, cpu_grid, Center(), Center(), Center())
                    levels[i, j, 1] = k
                    break
                end
            end
        end
    end

    return levels
end

@inline bottom_level(bbl, i, j) = @inbounds round(Int, bbl.bottom_index[i, j, 1])

@inline function bottom_tracer(c, bbl, i, j)
    k = max(bottom_level(bbl, i, j), 1)
    return @inbounds c[i, j, k]
end

# Transport across the face joining `(i, j)` to its neighbour, zero unless the shallower of the two bottom
# cells holds the denser water. NEMO's `zgdrho * mgrhu > 0`, where `zgdrho` is a buoyancy gradient
# (-δρ/ρ ≈ α δT - β δS) and `mgrhu` the sign of the bathymetry gradient. Only the sign of the product is
# used, so the α/β linearisation costs no accuracy against a full equation-of-state comparison.
@inline function active_face_transport(bbl, fields, grid, i, j, iⁿ, jⁿ, face_width, gradient_length)

    k  = bottom_level(bbl, i,  j)
    kⁿ = bottom_level(bbl, iⁿ, jⁿ)
    wet = (k > 0) & (kⁿ > 0)

    k⁺  = max(k,  1)
    kⁿ⁺ = max(kⁿ, 1)

    @inbounds begin
        T  = fields.T[i,  j,  k⁺];   S  = fields.S[i,  j,  k⁺]
        Tⁿ = fields.T[iⁿ, jⁿ, kⁿ⁺];  Sⁿ = fields.S[iⁿ, jⁿ, kⁿ⁺]
    end

    z  = znode(i,  j,  k⁺,  grid, Center(), Center(), Center())
    zⁿ = znode(iⁿ, jⁿ, kⁿ⁺, grid, Center(), Center(), Center())

    ℰ = bbl.equation_of_state
    α  = thermal_expansion(T,  S,  z,  ℰ);  β  = haline_contraction(T,  S,  z,  ℰ)
    αⁿ = thermal_expansion(Tⁿ, Sⁿ, zⁿ, ℰ);  βⁿ = haline_contraction(Tⁿ, Sⁿ, zⁿ, ℰ)

    δb = (αⁿ + α) * (Tⁿ - T) - (βⁿ + β) * (Sⁿ - S)
    slope = sign(z - zⁿ)                      # +1 where the bottom deepens toward the neighbour

    thickness = min(Δzᶜᶜᶜ(i, j, k⁺, grid), Δzᶜᶜᶜ(iⁿ, jⁿ, kⁿ⁺, grid))
    Q = bbl.diffusivity * face_width / gradient_length * thickness

    return ifelse(wet & (δb * slope > 0), Q, zero(grid))
end

@inline function x_face_transport(bbl::BottomBoundaryLayer, fields, grid, i, j)
    k = max(bottom_level(bbl, i, j), 1)
    return active_face_transport(bbl, fields, grid, i, j, i+1, j,
                                 Δyᶠᶜᶜ(i+1, j, k, grid), Δxᶠᶜᶜ(i+1, j, k, grid))
end

@inline function y_face_transport(bbl::BottomBoundaryLayer, fields, grid, i, j)
    k = max(bottom_level(bbl, i, j), 1)
    return active_face_transport(bbl, fields, grid, i, j, i, j+1,
                                 Δxᶜᶠᶜ(i, j+1, k, grid), Δyᶜᶠᶜ(i, j+1, k, grid))
end

@inline tracer_field(fields, ::Val{:T}) = fields.T
@inline tracer_field(fields, ::Val{:S}) = fields.S

"""
    bottom_boundary_layer_tendency(i, j, k, grid, clock, fields, parameters)

Discrete-form forcing. Nonzero only in the bottom cell; every other cell returns zero through `ifelse`.

Divides by `Az * Δz` where NEMO divides by area alone, because NEMO's tracer tendency is thickness
weighted and Oceananigans' is not. `Δz` is evaluated live so the metric follows z-star.
"""
@inline function bottom_boundary_layer_tendency(i, j, k, grid, clock, fields, parameters)

    bbl = parameters.bottom_boundary_layer
    c   = tracer_field(fields, parameters.tracer_name)

    kᵇ = bottom_level(bbl, i, j)
    on_bottom = (kᵇ > 0) & (k == kᵇ)
    kᵇ⁺ = max(kᵇ, 1)

    cᶜ = bottom_tracer(c, bbl, i, j)

    @inbounds begin
        Qe = bbl.transport_x[i,   j,   1]
        Qw = bbl.transport_x[i-1, j,   1]
        Qn = bbl.transport_y[i,   j,   1]
        Qs = bbl.transport_y[i,   j-1, 1]
    end

    divergence = Qe * (bottom_tracer(c, bbl, i+1, j) - cᶜ) +
                 Qw * (bottom_tracer(c, bbl, i-1, j) - cᶜ) +
                 Qn * (bottom_tracer(c, bbl, i, j+1) - cᶜ) +
                 Qs * (bottom_tracer(c, bbl, i, j-1) - cᶜ)

    volume = Azᶜᶜᶜ(i, j, kᵇ⁺, grid) * Δzᶜᶜᶜ(i, j, kᵇ⁺, grid)

    return ifelse(on_bottom, divergence / volume, zero(grid))
end

"""
    bottom_boundary_layer_forcing(bbl::BottomBoundaryLayer)

Return `(T, S)` forcings to merge into the ocean model's `forcing` keyword.
"""
function bottom_boundary_layer_forcing(bbl::BottomBoundaryLayer)
    T = Forcing(bottom_boundary_layer_tendency; discrete_form = true,
                parameters = (bottom_boundary_layer = bbl, tracer_name = Val(:T)))

    S = Forcing(bottom_boundary_layer_tendency; discrete_form = true,
                parameters = (bottom_boundary_layer = bbl, tracer_name = Val(:S)))

    return (; T, S)
end

"""
    BottomBoundaryLayerUpdate(bbl)

Callable that refreshes the face transports; register once per step with `add_callback!`.
"""
struct BottomBoundaryLayerUpdate{B}
    bottom_boundary_layer :: B
end

(u::BottomBoundaryLayerUpdate)(sim) = update_bottom_boundary_layer!(sim, u.bottom_boundary_layer)

"""
    bottom_boundary_layer_forcing(grid, diffusivity)

Convenience form for `omip_simulation`. Returns `(forcing, bbl)`; the forcing is an empty
`NamedTuple` and `bbl` is `nothing` when `diffusivity` is not set, so the scheme costs nothing when
it is off. The caller must register `BottomBoundaryLayerUpdate(bbl)` as a callback.
"""
function bottom_boundary_layer_forcing(grid, diffusivity)
    isnothing(diffusivity) && return NamedTuple(), nothing
    bbl = BottomBoundaryLayer(grid, TEOS10EquationOfState(); diffusivity)
    return bottom_boundary_layer_forcing(bbl), bbl
end
