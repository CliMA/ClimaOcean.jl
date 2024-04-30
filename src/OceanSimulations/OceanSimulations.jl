module OceanSimulations

<<<<<<< HEAD
export load_balanced_regional_grid, ocean_simulation

using Oceananigans.Units
using Oceananigans.Advection: TracerAdvection
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, inactive_node

=======
>>>>>>> origin/main
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
<<<<<<< HEAD
using Oceananigans.Operators

include("load_balanced_regional_grid.jl")

# Some defaults
default_free_surface(grid) = SplitExplicitFreeSurface(grid; cfl=0.7)
=======

# Some defualts
default_free_surface(grid) = SplitExplicitFreeSurface(cfl=0.7; grid)
>>>>>>> origin/main

function default_ocean_closure()
    mixing_length = MixingLength(Cᵇ=0.01)
    turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(Cᵂϵ=1.0)
    return CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)
end

<<<<<<< HEAD
default_momentum_advection() = VectorInvariant(; vorticity_scheme = WENO(; order = 9),
                                                  vertical_scheme = Centered(),
                                                divergence_scheme = WENO())

default_tracer_advection() = TracerAdvection(WENO(; order = 7),
                                             WENO(; order = 7),
                                             Centered())

@inline ϕ²(i, j, k, grid, ϕ)    = @inbounds ϕ[i, j, k]^2
@inline spᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline spᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * spᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_quadratic_bottom_drag(i, j, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * spᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_quadratic_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * spᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_quadratic_bottom_drag(i, j, k, grid, c, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * spᶜᶠᶜ(i, j, k, grid, Φ)

@inline is_immersed_drag_u(i, j, k, grid) = Int(immersed_peripheral_node(i, j, k-1, grid, Face(), Center(), Center()) & !inactive_node(i, j, k, grid, Face(), Center(), Center()))
@inline is_immersed_drag_v(i, j, k, grid) = Int(immersed_peripheral_node(i, j, k-1, grid, Center(), Face(), Center()) & !inactive_node(i, j, k, grid, Center(), Face(), Center()))

# Keep a constant linear drag parameter independent on vertical level
@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * is_immersed_drag_u(i, j, k, grid) * spᶠᶜᶜ(i, j, k, grid, fields) / Δzᶠᶜᶜ(i, j, k, grid)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * is_immersed_drag_v(i, j, k, grid) * spᶜᶠᶜ(i, j, k, grid, fields) / Δzᶜᶠᶜ(i, j, k, grid)

# TODO: Specify the grid to a grid on the sphere; otherwise we can provide a different
# function that requires latitude and longitude etc for computing coriolis=FPlane...
function ocean_simulation(grid; Δt = 5minutes,
=======
# TODO: Specify the grid to a grid on the sphere; otherwise we can provide a different
# function that requires latitude and longitude etc for computing coriolis=FPlane...
function ocean_simulation(grid;
>>>>>>> origin/main
                          closure = default_ocean_closure(),
                          free_surface = default_free_surface(grid),
                          reference_density = 1020,
                          rotation_rate = Ω_Earth,
<<<<<<< HEAD
                          gravitational_acceleration = g_Earth,
                          drag_coefficient = 0.003,
                          momentum_advection = default_momentum_advection(),
                          tracer_advection = default_tracer_advection(),
                          verbose = false)
=======
                          gravitational_acceleration = g_Earth)
>>>>>>> origin/main

    # Set up boundary conditions using Field
    top_zonal_momentum_flux      = Jᵘ = Field{Face, Center, Nothing}(grid)
    top_meridional_momentum_flux = Jᵛ = Field{Center, Face, Nothing}(grid)
    top_ocean_heat_flux          = Jᵀ = Field{Center, Center, Nothing}(grid)
    top_salt_flux                = Jˢ = Field{Center, Center, Nothing}(grid)

<<<<<<< HEAD
    u_bot_bc = FluxBoundaryCondition(u_quadratic_bottom_drag, discrete_form=true, parameters=drag_coefficient)
    v_bot_bc = FluxBoundaryCondition(v_quadratic_bottom_drag, discrete_form=true, parameters=drag_coefficient)

    ocean_boundary_conditions = (u = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵘ), bottom = u_bot_bc),
                                 v = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵛ), bottom = v_bot_bc),
                                 T = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ)),
                                 S = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ)))

    Fu = Forcing(u_immersed_bottom_drag, discrete_form=true, parameters=drag_coefficient)
    Fv = Forcing(v_immersed_bottom_drag, discrete_form=true, parameters=drag_coefficient)

    forcing = (; u = Fu, v = Fv)
=======
    ocean_boundary_conditions = (u = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵘ)),
                                 v = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵛ)),
                                 T = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵀ)),
                                 S = FieldBoundaryConditions(top=FluxBoundaryCondition(Jˢ)))
>>>>>>> origin/main

    # Use the TEOS10 equation of state
    teos10 = TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(; gravitational_acceleration, equation_of_state=teos10)

    # Minor simplifications for single column grids
<<<<<<< HEAD
    Nx, Ny, _ = size(grid)
    if Nx == Ny == 1 # single column grid
        tracer_advection = nothing
        momentum_advection = nothing
=======
    Nx, Ny, Nz = size(grid)
    if Nx == Ny == 1 # single column grid
        tracer_advection = nothing
        momentum_advection = nothing
    else
        # TODO: better advection scheme
        tracer_advection = WENO()
        momentum_advection = VectorInvariant(vorticity_scheme = WENO(),
                                             divergence_scheme = WENO(),
                                             vertical_scheme = WENO())
>>>>>>> origin/main
    end

    tracers = (:T, :S)
    if closure isa CATKEVerticalDiffusivity
        tracers = tuple(tracers..., :e)
<<<<<<< HEAD
        tracer_advection = (; T = tracer_advection, S = tracer_advection, e = nothing)
    end

    coriolis = HydrostaticSphericalCoriolis(; rotation_rate, scheme = ActiveCellEnstrophyConserving())
=======
    end

    coriolis = HydrostaticSphericalCoriolis(; rotation_rate)
>>>>>>> origin/main

    ocean_model = HydrostaticFreeSurfaceModel(; grid,
                                              buoyancy,
                                              closure,
                                              tracer_advection,
                                              momentum_advection,
                                              tracers,
                                              free_surface,
                                              coriolis,
<<<<<<< HEAD
                                              forcing,
                                              boundary_conditions = ocean_boundary_conditions)

    ocean = Simulation(ocean_model; Δt, verbose)
=======
                                              boundary_conditions = ocean_boundary_conditions)

    ocean = Simulation(ocean_model; Δt=5minutes, verbose=false)
>>>>>>> origin/main

    return ocean
end

end # module
