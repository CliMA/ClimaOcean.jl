using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities:
    CATKEVerticalDiffusivity,
    MixingLength,
    TurbulentKineticEnergyEquation

using Oceananigans.Grids: halo_size, topology, φnodes, λnodes
using Oceananigans.Fields: ConstantField, ZeroField

using ClimaSeaIce
using ClimaSeaIce.HeatBoundaryConditions: IceWaterThermalEquilibrium

using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

function regional_ecco2_grid(arch, Te, other_fields...; latitude, longitude)

    i₁ = 4 * first(longitude) + 1
    i₂ = 1440 - 4 * (360 - last(longitude))
    i₂ > i₁ || error("longitude $longitude is invalid.")
    Nx = i₂ - i₁ + 1

    j₁ = 4 * (90 + first(latitude)) + 1
    j₂ = 720 - 4 * (90 - last(latitude))
    j₂ > j₁ || error("latitude $latitude is invalid.")
    Ny = j₂ - j₁ + 1

    zc = znodes(Te)
    zf = znodes(Te.grid, Face())
    Δz = first(zspacings(Te.grid, Center()))

    Tᵢ = interior(Te, i₁:i₂, j₁:j₂, :)

    # Construct bottom_height depth by analyzing T
    Nx, Ny, Nz = size(Tᵢ)
    bottom_height = ones(Nx, Ny) .* (zf[1] - Δz)

    land = Tᵢ .< -10
    Tᵢ[land] .= NaN

    for i = 1:Nx, j = 1:Ny
        @inbounds for k = Nz:-1:1
            if isnan(Tᵢ[i, j, k])
                bottom_height[i, j] = zf[k+1]
                break
            end
        end
    end

    Tᵢ = arch_array(arch, Tᵢ)

    Tx = if longitude[2] - longitude[1] == 360
        Periodic
    else
        Bounded
    end

    grid = LatitudeLongitudeGrid(arch; latitude, longitude,
                                 size = (Nx, Ny, Nz),
                                 halo = (7, 7, 7),
                                 z = zf,
                                 topology = (Periodic, Bounded, Bounded))

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

    Nf = length(other_fields)

    ft = ntuple(Nf) do n
        fe = other_fields[n]
        fᵢ = interior(fe, i₁:i₂, j₁:j₂, :)
        fᵢ[land] .= NaN
        fᵢ = arch_array(arch, fᵢ)
    end

    all_fields = tuple(Tᵢ, ft...)

    return grid, all_fields
end

function omip_ocean_component(grid;
                              closure = :default,
                              passive_tracers = tuple())

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

    tracers = tuple(:T, :S, passive_tracers...)

    if closure == :default
        mixing_length = MixingLength(Cᵇ=0.01)
        turbulent_kinetic_energy_equation = TurbulentKineticEnergyEquation(Cᵂϵ=1.0)
        closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation)
        tracers = tuple(:e, tracers...)
    end
    
    ocean_model = HydrostaticFreeSurfaceModel(; grid, buoyancy, closure,
                                              tracer_advection, momentum_advection,
                                              tracers = (:T, :S), #, :e),
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
                                    advection = WENO(),
                                    ice_consolidation_thickness = 0.05,
                                    ice_salinity = 4,
                                    internal_heat_flux = ConductiveFlux(conductivity=2),
                                    top_heat_flux = ConstantField(0), # W m⁻²
                                    top_heat_boundary_condition = PrescribedTemperature(-10),
                                    bottom_heat_boundary_condition = bottom_bc,
                                    bottom_heat_flux = sea_ice_ocean_heat_flux)

    sea_ice = Simulation(sea_ice_model, Δt=5minutes, verbose=false)

    elapsed = time_ns() - start_time
    msg = string("Finished building sea ice component. (" * prettytime(elapsed * 1e-9), ")")
    @info msg

    return sea_ice
end

