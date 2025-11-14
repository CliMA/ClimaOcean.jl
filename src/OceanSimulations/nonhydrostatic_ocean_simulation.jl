const default_nonhydrostatic_domain = (Lx = 256meters, Lz = 128meters)
const default_nonhydrostatic_spacing = (Δx = 4meters, Δz = 2meters)
const default_nonhydrostatic_size = (
    Nx = Int(round(default_nonhydrostatic_domain.Lx / default_nonhydrostatic_spacing.Δx)),
    Nz = Int(round(default_nonhydrostatic_domain.Lz / default_nonhydrostatic_spacing.Δz)),
)

"""
    default_nonhydrostatic_grid(arch=CPU())

Return a 2D RectilinearGrid (x–z slice) with 4 m horizontal spacing over a 256 m
domain and 2 m vertical spacing over 128 m depth, configured for LES experiments.
"""
function default_nonhydrostatic_grid(arch=CPU())
    Nx = default_nonhydrostatic_size.Nx
    Nz = default_nonhydrostatic_size.Nz
    Ny = 1
    Ly = default_nonhydrostatic_spacing.Δx # retain a thin y-direction for 2D slice

    grid = RectilinearGrid(arch;
                           size = (Nx, Ny, Nz),
                           x = (0, default_nonhydrostatic_domain.Lx),
                           y = (0, Ly),
                           z = (-default_nonhydrostatic_domain.Lz, 0),
                           topology = (Periodic, Flat, Bounded))

    return grid
end

"""
    nonhydrostatic_ocean_simulation(; kwargs...)

Build a NonhydrostaticModel-based ocean simulation configured for high-resolution LES.
The default configuration follows a 2D nighttime boundary-layer experiment with WENO
advection, TEOS-10 equation of state, and no turbulence closure.
"""
function nonhydrostatic_ocean_simulation(grid; Δt=1,
                                         tracers = (:T, :S),
                                         gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration,
                                         coriolis = nothing,
                                         closure = nothing,
                                         forcing = NamedTuple(),
                                         biogeochemistry = nothing,
                                         advection = WENO(order=9),
                                         equation_of_state = TEOS10EquationOfState(),
                                         boundary_conditions::NamedTuple = NamedTuple(),
                                         verbose = false)

    FT = eltype(grid)

    τx = Field{Face, Center, Nothing}(grid)
    τy = Field{Center, Face, Nothing}(grid)
    Jᵀ = Field{Center, Center, Nothing}(grid)
    Jˢ = Field{Center, Center, Nothing}(grid)

    u_top_bc = FluxBoundaryCondition(τx)
    v_top_bc = FluxBoundaryCondition(τy)
    T_top_bc = FluxBoundaryCondition(Jᵀ)
    S_top_bc = FluxBoundaryCondition(Jˢ)

    default_boundary_conditions = (u = FieldBoundaryConditions(top=u_top_bc),
                                   v = FieldBoundaryConditions(top=v_top_bc),
                                   T = FieldBoundaryConditions(top=T_top_bc),
                                   S = FieldBoundaryConditions(top=S_top_bc))

    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    buoyancy = SeawaterBuoyancy(; gravitational_acceleration, equation_of_state)

    ocean_model = NonhydrostaticModel(; grid,
                                      buoyancy,
                                      closure,
                                      biogeochemistry,
                                      advection,
                                      tracers,
                                      coriolis,
                                      forcing,
                                      boundary_conditions)

    ocean = Simulation(ocean_model; Δt, verbose)
    conjure_time_step_wizard!(ocean, cfl=0.7)

    return ocean
end

