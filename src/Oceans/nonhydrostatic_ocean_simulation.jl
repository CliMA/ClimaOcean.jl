
"""
    nonhydrostatic_ocean_simulation(grid; kwargs...)

Build a NonhydrostaticModel-based ocean simulation configured for high-resolution LES.
The default configuration has WENO advection, TEOS-10 equation of state, and no turbulence closure.
"""
function nonhydrostatic_ocean_simulation(grid; Δt=1,
                                         tracers = (:T, :S),
                                         gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration,
                                         coriolis = nothing,
                                         stokes_drift = nothing,
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
                                      stokes_drift,
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

