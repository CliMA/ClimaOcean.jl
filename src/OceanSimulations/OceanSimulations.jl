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

include("load_balanced_regional_grid.jl")

# Some defualts
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

# TODO: Specify the grid to a grid on the sphere; otherwise we can provide a different
# function that requires latitude and longitude etc for computing coriolis=FPlane...
function ocean_simulation(grid; Δt = 5minutes,
                          closure = default_ocean_closure(),
                          free_surface = default_free_surface(grid),
                          reference_density = 1020,
                          rotation_rate = Ω_Earth,
                          gravitational_acceleration = g_Earth,
                          momentum_advection = default_momentum_advection(),
                          tracer_advection = default_tracer_advection())

    # Set up boundary conditions using Field
    top_zonal_momentum_flux      = Jᵘ = Field{Face, Center, Nothing}(grid)
    top_meridional_momentum_flux = Jᵛ = Field{Center, Face, Nothing}(grid)
    top_ocean_heat_flux          = Jᵀ = Field{Center, Center, Nothing}(grid)
    top_salt_flux                = Jˢ = Field{Center, Center, Nothing}(grid)

    ocean_boundary_conditions = (u = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵘ)),
                                 v = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵛ)),
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
