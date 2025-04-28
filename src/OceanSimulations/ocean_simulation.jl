using Oceananigans.DistributedComputations: DistributedGrid, all_reduce
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: DefaultBoundaryCondition
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node, MutableGridOfSomeKind
using Oceananigans.OrthogonalSphericalShellGrids

using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    CATKEMixingLength,
    CATKEEquation

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Statistics: mean

# Some defaults
default_free_surface(grid) = SplitExplicitFreeSurface(grid; cfl=0.7)

estimate_maximum_Δt(grid::RectilinearGrid) = 30minutes # ?

function estimate_maximum_Δt(grid)
    arch = architecture(grid)
    Δx = mean(xspacings(grid))
    Δy = mean(yspacings(grid))
    Δθ = rad2deg(mean([Δx, Δy])) / grid.radius

    # The maximum Δt is roughly 30minutes / Δθ, giving:
    # - 30 minutes for a 1 degree ocean
    # - 15 minutes for a 1/4 degree ocean
    # - 7.5 minutes for a 1/8 degree ocean
    # - 3.75 minutes for a 1/16 degree ocean
    # - 1.875 minutes for a 1/32 degree ocean

    Δt = 30minutes / Δθ

    return all_reduce(min, Δt, arch)
end

const TripolarOfSomeKind = Union{TripolarGrid, ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:TripolarGrid}}

function default_free_surface(grid::TripolarOfSomeKind;
                              fixed_Δt = estimate_maximum_Δt(grid),
                              cfl = 0.7)
    free_surface = SplitExplicitFreeSurface(grid; cfl, fixed_Δt)
    return free_surface
end

function default_free_surface(grid::DistributedGrid;
                              fixed_Δt = estimate_maximum_Δt(grid),
                              cfl = 0.7)

    free_surface = SplitExplicitFreeSurface(grid; cfl, fixed_Δt)
    substeps = length(free_surface.substepping.averaging_weights)
    substeps = all_reduce(max, substeps, architecture(grid))
    free_surface = SplitExplicitFreeSurface(grid; substeps)
    @info "Using a $(free_surface)"
    return free_surface
end

default_vertical_coordinate(grid) = Oceananigans.Models.ZCoordinate()
default_vertical_coordinate(::MutableGridOfSomeKind) = Oceananigans.Models.ZStar()

function default_ocean_closure(FT=Oceananigans.defaults.FloatType)
    mixing_length = CATKEMixingLength(Cᵇ=0.01)
    turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ=1.0)

    return CATKEVerticalDiffusivity(VerticallyImplicitTimeDiscretization(), FT; mixing_length, turbulent_kinetic_energy_equation)
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

# Keep a constant linear drag parameter independent on vertical level
@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, fields)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, fields)

# TODO: Specify the grid to a grid on the sphere; otherwise we can provide a different
# function that requires latitude and longitude etc for computing coriolis=FPlane...
function ocean_simulation(grid;
                          Δt = estimate_maximum_Δt(grid),
                          closure = default_ocean_closure(),
                          tracers = (:T, :S),
                          free_surface = default_free_surface(grid),
                          reference_density = 1020,
                          rotation_rate = Ω_Earth,
                          gravitational_acceleration = g_Earth,
                          bottom_drag_coefficient = Default(0.003),
                          forcing = NamedTuple(),
                          biogeochemistry = nothing,
                          timestepper = :QuasiAdamsBashforth2,
                          coriolis = Default(HydrostaticSphericalCoriolis(; rotation_rate)),
                          momentum_advection = default_momentum_advection(),
                          equation_of_state = TEOS10EquationOfState(; reference_density),
                          boundary_conditions::NamedTuple = NamedTuple(),
                          tracer_advection = default_tracer_advection(),
                          vertical_coordinate = default_vertical_coordinate(grid),
                          warn = true,
                          verbose = false)

    FT = eltype(grid)

    if grid isa RectilinearGrid # turn off Coriolis unless user-supplied
        coriolis = default_or_override(coriolis, nothing)
    else
        coriolis = default_or_override(coriolis)
    end

    # Detect whether we are on a single column grid
    Nx, Ny, _ = size(grid)
    single_column_simulation = Nx == 1 && Ny == 1

    if single_column_simulation
        # Let users put a bottom drag if they want
        bottom_drag_coefficient = default_or_override(bottom_drag_coefficient, zero(grid))

        # Don't let users use advection in a single column model
        tracer_advection = nothing
        momentum_advection = nothing

        # No immersed boundaries in a single column grid
        u_immersed_bc = DefaultBoundaryCondition()
        v_immersed_bc = DefaultBoundaryCondition()
    else
        if warn && !(grid isa ImmersedBoundaryGrid)
            msg = """Are you totally, 100% sure that you want to build a simulation on

                   $(summary(grid))

                   rather than on an ImmersedBoundaryGrid?
                   """
            @warn msg
        end

        bottom_drag_coefficient = default_or_override(bottom_drag_coefficient)

        u_immersed_drag = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
        v_immersed_drag = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)

        u_immersed_bc = ImmersedBoundaryCondition(bottom=u_immersed_drag)
        v_immersed_bc = ImmersedBoundaryCondition(bottom=v_immersed_drag)

        # Forcing for u, v
        u_barotropic_potential = Field{Center, Center, Nothing}(grid)
        v_barotropic_potential = Field{Center, Center, Nothing}(grid)
        u_forcing = BarotropicPotentialForcing(XDirection(), u_barotropic_potential)
        v_forcing = BarotropicPotentialForcing(YDirection(), v_barotropic_potential)

        :u ∈ keys(forcing) && (u_forcing = (u_forcing, forcing[:u]))
        :v ∈ keys(forcing) && (v_forcing = (v_forcing, forcing[:v]))
        forcing = merge(forcing, (u=u_forcing, v=v_forcing))
    end

    bottom_drag_coefficient = convert(FT, bottom_drag_coefficient)

    # Set up boundary conditions using Field
    top_zonal_momentum_flux      = τx = Field{Face, Center, Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center, Face, Nothing}(grid)
    top_ocean_heat_flux          = Jᵀ = Field{Center, Center, Nothing}(grid)
    top_salt_flux                = Jˢ = Field{Center, Center, Nothing}(grid)

    # Construct ocean boundary conditions including surface forcing and bottom drag
    u_top_bc = FluxBoundaryCondition(τx)
    v_top_bc = FluxBoundaryCondition(τy)
    T_top_bc = FluxBoundaryCondition(Jᵀ)
    S_top_bc = FluxBoundaryCondition(Jˢ)

    u_bot_bc = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)
    v_bot_bc = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form=true, parameters=bottom_drag_coefficient)

    default_boundary_conditions = (u = FieldBoundaryConditions(top=u_top_bc, bottom=u_bot_bc, immersed=u_immersed_bc),
                                   v = FieldBoundaryConditions(top=v_top_bc, bottom=v_bot_bc, immersed=v_immersed_bc),
                                   T = FieldBoundaryConditions(top=T_top_bc),
                                   S = FieldBoundaryConditions(top=S_top_bc))

    # Merge boundary conditions with preference to user
    # TODO: support users specifying only _part_ of the bcs for u, v, T, S (ie adding the top and immersed
    # conditions even when a user-bc is supplied).
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
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
                                              biogeochemistry,
                                              tracer_advection,
                                              momentum_advection,
                                              tracers,
                                              timestepper,
                                              free_surface,
                                              coriolis,
                                              forcing,
                                              boundary_conditions,
                                              vertical_coordinate)

    ocean = Simulation(ocean_model; Δt, verbose)

    return ocean
end

hasclosure(closure, ClosureType) = closure isa ClosureType
hasclosure(closure_tuple::Tuple, ClosureType) = any(hasclosure(c, ClosureType) for c in closure_tuple)

