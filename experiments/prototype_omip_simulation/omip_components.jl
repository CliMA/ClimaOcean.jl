using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation

using Oceananigans.Grids: halo_size, topology, φnodes, λnodes
using Oceananigans.Fields: ConstantField, ZeroField

using ClimaSeaIce
using ClimaSeaIce.HeatBoundaryConditions: IceWaterThermalEquilibrium

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

function omip_ocean_component(grid)

    start_time = time_ns()

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

    Nx, Ny, Nz = size(grid)

    if Nx == Ny == 1
        tracer_advection = nothing
        momentum_advection = nothing
    else
        tracer_advection = WENO()
        momentum_advection = VectorInvariant(vorticity_scheme = WENO(),
                                             divergence_scheme = WENO(),
                                             vertical_scheme = WENO())
    end

    mixing_length = MixingLength(Cᵇ=0.01)
    turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(Cᵂϵ=1.0)
    closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)
    
    ocean_model = HydrostaticFreeSurfaceModel(; grid, buoyancy, closure,
                                              tracer_advection, momentum_advection,
                                              tracers = (:T, :S, :e),
                                              free_surface = SplitExplicitFreeSurface(cfl=0.7; grid),
                                              boundary_conditions = ocean_boundary_conditions,
                                              coriolis = HydrostaticSphericalCoriolis())

    ocean = Simulation(ocean_model; Δt=5minutes, verbose=false)

    elapsed = time_ns() - start_time
    msg = string("Finished building ocean component. (" * prettytime(elapsed * 1e-9), ")")
    @info msg

    return ocean
end

function omip_sea_ice_component(ocean_model)
    start_time = time_ns()

    ocean_grid = ocean_model.grid
    Nx, Ny, Nz = size(ocean_grid)
    Hx, Hy, Hz = halo_size(ocean_grid)
    TX, TY, TZ = topology(ocean_grid)

    λ = λnodes(ocean_grid, Face())
    φ = φnodes(ocean_grid, Face())
    longitude = (λ[1], λ[end])
    latitude  = (φ[1], φ[end])

    sea_ice_grid = LatitudeLongitudeGrid(arch; longitude, latitude,
                                         size = (Nx, Ny),
                                         halo = (Hx, Hy),
                                         topology = (TX, TY, Flat))

    if ocean_grid isa ImmersedBoundaryGrid
        h = ocean_grid.immersed_boundary.bottom_height
        land = interior(h) .>= 0
        sea_ice_grid = ImmersedBoundaryGrid(sea_ice_grid, GridFittedBoundary(land))
    end

    sea_ice_ocean_heat_flux = Field{Center, Center, Nothing}(ocean_grid)

    Nz = size(ocean_grid, 3)
    So = ocean_model.tracers.S
    ocean_surface_salinity = view(So, :, :, Nz)
    bottom_bc = IceWaterThermalEquilibrium(ocean_surface_salinity)

    u, v, w = ocean_model.velocities
    ocean_surface_velocities = (u = view(u, :, :, Nz), #interior(u, :, :, Nz),
                                v = view(v, :, :, Nz), #interior(v, :, :, Nz),
                                w = ZeroField())

    sea_ice_model = SlabSeaIceModel(sea_ice_grid;
                                    velocities = ocean_surface_velocities,
                                    advection = nothing,
                                    ice_consolidation_thickness = 0.05,
                                    ice_salinity = 4,
                                    internal_heat_flux = ConductiveFlux(conductivity=2),
                                    top_heat_flux = ConstantField(0), # W m⁻²
                                    top_heat_boundary_condition = PrescribedTemperature(0),
                                    bottom_heat_boundary_condition = bottom_bc,
                                    bottom_heat_flux = sea_ice_ocean_heat_flux)

    sea_ice = Simulation(sea_ice_model, Δt=5minutes, verbose=false)

    elapsed = time_ns() - start_time
    msg = string("Finished building sea ice component. (" * prettytime(elapsed * 1e-9), ")")
    @info msg

    return sea_ice
end

