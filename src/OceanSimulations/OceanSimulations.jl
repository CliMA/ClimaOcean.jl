module OceanSimulations

export ocean_simulation

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: with_tracers
using Oceananigans.Advection: FluxFormAdvection
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node
using OrthogonalSphericalShellGrids

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    CATKEMixingLength,
    CATKEEquation

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Operators

struct Default{V}
    value :: V
end

Default() = Default(nothing)
default_or_override(default::Default, value=default.value) = value

# Some defaults
default_free_surface(grid) = SplitExplicitFreeSurface(grid; cfl=0.7)

# 70 substeps is a safe rule of thumb for an ocean at 1/4 - 1/10th of a degree
# TODO: pass the cfl and a given Δt to calculate the number of substeps?
default_free_surface(grid::TripolarGrid) = SplitExplicitFreeSurface(grid; substeps = 70)

function default_ocean_closure()
    mixing_length = CATKEMixingLength(Cᵇ=0.01)
    turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ=1.0)
    return CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)
end

default_momentum_advection() = VectorInvariant(; vorticity_scheme = WENO(order=9),
                                                  vertical_scheme = Centered(),
                                                divergence_scheme = WENO(order=5))

default_tracer_advection() = FluxFormAdvection(WENO(order=7),
                                               WENO(order=7),
                                               Centered())

@inline ϕ²(i, j, k, grid, ϕ)    = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline is_immersed_drag_u(i, j, k, grid) = Int(immersed_peripheral_node(i, j, k-1, grid, Face(), Center(), Center()) & !inactive_node(i, j, k, grid, Face(), Center(), Center()))
@inline is_immersed_drag_v(i, j, k, grid) = Int(immersed_peripheral_node(i, j, k-1, grid, Center(), Face(), Center()) & !inactive_node(i, j, k, grid, Center(), Face(), Center()))

# Keep a constant linear drag parameter independent on vertical level
@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * is_immersed_drag_u(i, j, k, grid) * spᶠᶜᶜ(i, j, k, grid, fields) / Δzᶠᶜᶜ(i, j, k, grid)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * is_immersed_drag_v(i, j, k, grid) * spᶜᶠᶜ(i, j, k, grid, fields) / Δzᶜᶠᶜ(i, j, k, grid)

# TODO: Specify the grid to a grid on the sphere; otherwise we can provide a different
# function that requires latitude and longitude etc for computing coriolis=FPlane...
function ocean_simulation(grid;
                          Δt = 5minutes,
                          closure = default_ocean_closure(),
                          tracers = (:T, :S),
                          free_surface = default_free_surface(grid),
                          reference_density = 1020,
                          rotation_rate = Ω_Earth,
                          gravitational_acceleration = g_Earth,
                          bottom_drag_coefficient = Default(0.003),
                          forcing = NamedTuple(),
                          coriolis = HydrostaticSphericalCoriolis(; rotation_rate),
                          momentum_advection = default_momentum_advection(),
                          equation_of_state = TEOS10EquationOfState(; reference_density),
                          tracer_advection = default_tracer_advection(),
                          verbose = false)

    FT = eltype(grid)

    # Detect whether we are on a single column grid
    Nx, Ny, _ = size(grid)
    single_column_simulation = Nx == 1 && Ny == 1

    if single_column_simulation
        # Let users put a bottom drag if they want
        bottom_drag_coefficient = default_or_override(bottom_drag_coefficient, zero(grid))

        # Don't let users use advection in a single column model
        tracer_advection = nothing
        momentum_advection = nothing
    else
        bottom_drag_coefficient = default_or_override(bottom_drag_coefficient)
    end

    # Set up boundary conditions using Field
    top_zonal_momentum_flux      = τx = Field{Face, Center, Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center, Face, Nothing}(grid)
    top_ocean_heat_flux          = Jᵀ = Field{Center, Center, Nothing}(grid)
    top_salt_flux                = Jˢ = Field{Center, Center, Nothing}(grid)

    u_bot_bc = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    v_bot_bc = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)

    ocean_boundary_conditions = (u = FieldBoundaryConditions(top = FluxBoundaryCondition(τx), bottom = u_bot_bc),
                                 v = FieldBoundaryConditions(top = FluxBoundaryCondition(τy), bottom = v_bot_bc),
                                 T = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ)),
                                 S = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ)))

    if grid isa ImmersedBoundaryGrid
        Fu = Forcing(u_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
        Fv = Forcing(v_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
        forcing = merge(forcing, (u=Fu, v=Fv))
    end

    buoyancy = SeawaterBuoyancy(; gravitational_acceleration, equation_of_state)

    if tracer_advection isa NamedTuple
        tracer_advection = with_tracers(tracers, tracer_advection, default_tracer_advection())
    else
        tracer_advection = NamedTuple(name => tracer_advection for name in tracers)
    end

    if hasclosure(closure, CATKEVerticalDiffusivity)
        # Magically add :e to tracers
        if !(:e ∈ tracers)
            tracers = tuple(tracers..., :e)
        end

        # Turn off CATKE tracer advection
        tke_advection = (; e=nothing)
        tracer_advection = merge(tracer_advection, tke_advection)
    end

    ocean_model = HydrostaticFreeSurfaceModel(; grid,
                                              buoyancy,
                                              closure,
                                              tracer_advection,
                                              momentum_advection,
                                              tracers,
                                              free_surface,
                                              coriolis,
                                              forcing,
                                              boundary_conditions = ocean_boundary_conditions)

    ocean = Simulation(ocean_model; Δt, verbose)

    return ocean
end

hasclosure(closure, ClosureType) = closure isa ClosureType
hasclosure(closure_tuple::Tuple, ClosureType) = any(hasclosure(c, ClosureType) for c in closure_tuple)

end # module
