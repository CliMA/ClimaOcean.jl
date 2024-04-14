module OceanSimulations

export load_balanced_regional_grid, ocean_simulation

using Oceananigans.Units
using Oceananigans.Advection: TracerAdvection

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.Operators

include("load_balanced_regional_grid.jl")

# Some defaults
default_free_surface(grid) = SplitExplicitFreeSurface(; cfl=0.7, grid)

function default_ocean_closure()
    mixing_length = MixingLength(Cᵇ=0.01)
    turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(Cᵂϵ=1.0)
    return CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)
end

default_momentum_advection() = VectorInvariant(; vorticity_scheme = WENO(; order = 9),
                                                  vertical_scheme = Centered(),
                                                divergence_scheme = WENO())

default_tracer_advection() = TracerAdvection(WENO(; order = 7),
                                             WENO(; order = 7),
                                             Centered())

@inline ϕ²(i, j, k, grid, ϕ)       = @inbounds ϕ[i, j, k]^2
@inline speedᶠᶜᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², Φ.v))
@inline speedᶜᶠᶜ(i, j, k, grid, Φ) = @inbounds sqrt(Φ.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², Φ.u))

@inline u_drag_bc(i, j, grid, clock, Φ, μ) = @inbounds - μ * Φ.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, Φ)
@inline v_drag_bc(i, j, grid, clock, Φ, μ) = @inbounds - μ * Φ.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, Φ)

@inline u_immersed_drag_bc(i, j, k, grid, clock, Φ, μ) = @inbounds - μ * Φ.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, Φ)
@inline v_immersed_drag_bc(i, j, k, grid, clock, Φ, μ) = @inbounds - μ * Φ.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, Φ)

# TODO: Specify the grid to a grid on the sphere; otherwise we can provide a different
# function that requires latitude and longitude etc for computing coriolis=FPlane...
function ocean_simulation(grid; Δt = 5minutes,
                          closure = default_ocean_closure(),
                          free_surface = default_free_surface(grid),
                          reference_density = 1020,
                          rotation_rate = Ω_Earth,
                          gravitational_acceleration = g_Earth,
                          drag_coefficient = 0.003,
                          momentum_advection = default_momentum_advection(),
                          tracer_advection = default_tracer_advection())

    # Set up boundary conditions using Field
    top_zonal_momentum_flux      = Jᵘ = Field{Face, Center, Nothing}(grid)
    top_meridional_momentum_flux = Jᵛ = Field{Center, Face, Nothing}(grid)
    top_ocean_heat_flux          = Jᵀ = Field{Center, Center, Nothing}(grid)
    top_salt_flux                = Jˢ = Field{Center, Center, Nothing}(grid)

    u_bottom_drag   = FluxBoundaryCondition(u_drag_bc, discrete_form=true, parameters=drag_coefficient)
    v_bottom_drag   = FluxBoundaryCondition(v_drag_bc, discrete_form=true, parameters=drag_coefficient)
    u_immersed_drag = FluxBoundaryCondition(u_immersed_drag_bc, discrete_form=true, parameters=drag_coefficient)
    v_immersed_drag = FluxBoundaryCondition(v_immersed_drag_bc, discrete_form=true, parameters=drag_coefficient)

    u_immersed_bc = ImmersedBoundaryCondition(bottom=u_immersed_drag)
    v_immersed_bc = ImmersedBoundaryCondition(bottom=v_immersed_drag)

    ocean_boundary_conditions = (u = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵘ), bottom = u_bottom_drag, immersed = u_immersed_bc),
                                 v = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵛ), bottom = v_bottom_drag, immersed = v_immersed_bc),
                                 T = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵀ)),
                                 S = FieldBoundaryConditions(top=FluxBoundaryCondition(Jˢ)))

    # Use the TEOS10 equation of state
    teos10 = TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(; gravitational_acceleration, equation_of_state=teos10)

    # Minor simplifications for single column grids
    Nx, Ny, _ = size(grid)
    if Nx == Ny == 1 # single column grid
        tracer_advection = nothing
        momentum_advection = nothing
    end

    tracers = (:T, :S)
    if closure isa CATKEVerticalDiffusivity
        tracers = tuple(tracers..., :e)
        tracer_advection = (; T = tracer_advection, S = tracer_advection, e = nothing)
    end

    coriolis = HydrostaticSphericalCoriolis(; rotation_rate)

    ocean_model = HydrostaticFreeSurfaceModel(; grid,
                                              buoyancy,
                                              closure,
                                              tracer_advection,
                                              momentum_advection,
                                              tracers,
                                              free_surface,
                                              coriolis,
                                              boundary_conditions = ocean_boundary_conditions)

    ocean = Simulation(ocean_model; Δt, verbose=false)

    return ocean
end

end # module
