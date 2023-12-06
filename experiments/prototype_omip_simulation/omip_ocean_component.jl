top_ocean_heat_flux          = Qᵀ = Field{Center, Center, Nothing}(grid)
top_salt_flux                = Fˢ = Field{Center, Center, Nothing}(grid)
top_zonal_momentum_flux      = τˣ = Field{Face, Center, Nothing}(grid)
top_meridional_momentum_flux = τʸ = Field{Center, Face, Nothing}(grid)

ocean_boundary_conditions = (u = FieldBoundaryConditions(top=FluxBoundaryCondition(τˣ)),
                             v = FieldBoundaryConditions(top=FluxBoundaryCondition(τʸ)),
                             T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                             S = FieldBoundaryConditions(top=FluxBoundaryCondition(Fˢ)))

# Model construction
teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)

using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: MixingLength
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: TurbulentKineticEnergyEquation

mixing_length = MixingLength(Cᵇ=0.01)
turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(Cᵂϵ=1.0)
closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)

tracer_advection = WENO()
momentum_advection = VectorInvariant(vorticity_scheme = WENO(),
                                     divergence_scheme = WENO(),
                                     vertical_scheme = WENO())

ocean_model = HydrostaticFreeSurfaceModel(; grid, buoyancy, closure,
                                          tracer_advection, momentum_advection,
                                          tracers = (:T, :S, :e),
                                          free_surface = SplitExplicitFreeSurface(cfl=0.7; grid),
                                          boundary_conditions = ocean_boundary_conditions,
                                          coriolis = HydrostaticSphericalCoriolis())

ocean = Simulation(ocean_model; Δt=5minutes, verbose=false)

