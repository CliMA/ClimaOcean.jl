# Neverworld
#
# Ingredients:
#
#   * Zonally-periodic domain with continental shelves on all boundaries except Southern Ocean
#       * longitude = (0, 60)
#       * 2 configurations in latitude
#           - Half basin: latitude = (-70, 0)
#           - Full basin: latitude = (-70, 70)
#       * z to -4000m
#       * Southern Ocean channel from -60 to -40 (with a ridge to -2000 m)
#
#   * Zonally-homogeneous wind stress with mid-latitude jet and trade winds
#
#   * Buoyancy 
#       * restoring hot at the equator and cold at the poles (parabola, cosine, smooth step function)
#       * equator-pole buoyancy differential: 0.06 (α * g * 30 ≈ 0.06 with α=2e-4, g=9.81)
#       * exponential initial vertical stratification with N² = 6e-5 and decay scale h = 1000 m
#           - eg bᵢ = N² * h * exp(z / h)
#
#   * Quadratic bottom drag with drag_coefficient = 1e-3

using CubicSplines

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: ynode
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using ClimaOcean: CubicSplineFunction, u_bottom_drag, v_bottom_drag, u_immersed_bottom_drag, v_immersed_bottom_drag
using ..VerticalGrids: stretched_vertical_cell_interfaces

default_z = stretched_vertical_cell_interfaces(surface_layer_Δz = 5,
                                               surface_layer_height = 100,
                                               stretching_exponent = 1.02,
                                               minimum_depth = 3000)

# Default Neverworld wind stress profile
zonal_stress_latitudes = [-70,   -45,   -15,      0,    15,    45,  70]
zonal_stress_nodes     = [ +0,  -0.2,  +0.1,  +0.02,  +0.1,  -0.1,  +0] .* 1e-3 # kinematic wind stress
default_zonal_wind_stress = CubicSplineFunction{:y}(zonal_stress_latitudes, zonal_stress_nodes)

const c = Center()

@inline default_target_buoyancy_distribution(φ, p) = p.Δb * cos(π * φ / p.Δφ)

@inline function buoyancy_relaxation(i, j, grid, clock, fields, parameters)
    k = grid.Nz

    # Target buoyancy distribution
    φ = ynode(c, c, c, i, j, k, grid)#, c, c, c)
    b★ = parameters.b★(φ, parameters)

    Δz = Δzᶜᶜᶜ(i, j, k, grid)
    t★ = parameters.t★
    q★ = Δz / t★

    return @inbounds q★ * (fields.b[i, j, k] - b★)
end

function neverworld_simulation(arch;
                               size = (70, 60, length(default_z)-1), 
                               latitude = (-70, 0),
                               longitude = (0, 60),
                               z = default_z,
                               momentum_advection = VectorInvariant(vorticity_scheme   = WENO(),
                                                                    divergence_scheme  = WENO(),
                                                                    vertical_scheme    = WENO()),
                               tracer_advection = WENO(),
                               closure = CATKEVerticalDiffusivity(),
                               tracers = (:b, :e),
                               buoyancy = BuoyancyTracer(),
                               buoyancy_relaxation_time_scale = 4days,
                               target_buoyancy_distribution = default_target_buoyancy_distribution,
                               bottom_drag_coefficient = 2e-3,
                               equator_pole_buoyancy_difference = 0.06,
                               surface_buoyancy_gradient = 1e-5,
                               stratification_scale_height = 1000, # meters
                               time_step = 5minutes,
                               zonal_wind_stress = default_zonal_wind_stress)

    grid = LatitudeLongitudeGrid(arch; size, latitude, longitude, z, halo=(5, 5, 5))

    N² = surface_buoyancy_gradient
    h  = stratification_scale_height 
    Δb = equator_pole_buoyancy_difference 
    t★ = buoyancy_relaxation_time_scale 
    Δφ = - latitude[1]
    b★ = target_buoyancy_distribution 
    μ  = bottom_drag_coefficient

    # Buoyancy flux
    parameters = (; Δφ, Δb, t★, b★)
    b_top_bc = FluxBoundaryCondition(buoyancy_relaxation, discrete_form=true; parameters)

    # Wind stress
    zonal_wind_stress_field = Field{Face, Center, Nothing}(grid)
    set!(zonal_wind_stress_field, zonal_wind_stress) 
    u_top_bc = FluxBoundaryCondition(interior(zonal_wind_stress_field, :, :, 1))

    # Bottom drag
    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters=μ)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters=μ)

    u_immersed_bc = ImmersedBoundaryCondition(bottom=drag_u) 
    v_immersed_bc = ImmersedBoundaryCondition(bottom=drag_v) 

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form=true, parameters=μ)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form=true, parameters=μ)

    u_bcs = FieldBoundaryConditions(bottom=u_bottom_drag_bc, immersed=u_immersed_bc, top=u_top_bc)
    v_bcs = FieldBoundaryConditions(bottom=v_bottom_drag_bc, immersed=v_immersed_bc)
    b_bcs = FieldBoundaryConditions(top=b_top_bc)

    coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConservingScheme())
    free_surface = ImplicitFreeSurface()
    model = HydrostaticFreeSurfaceModel(; grid, tracers, buoyancy, coriolis, free_surface,
                                        momentum_advection, tracer_advection,
                                        boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs))

    bᵢ(x, y, z) = Δb + N² * h * (exp(z / h) - 1)
    set!(model, b=bᵢ)

    simulation = Simulation(model; Δt=time_step)

    return simulation
end

