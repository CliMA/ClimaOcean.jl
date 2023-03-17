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
using Oceananigans.Grids: xnode, ynode
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity

using ClimaOcean: CubicSplineFunction, u_bottom_drag, v_bottom_drag, u_immersed_bottom_drag, v_immersed_bottom_drag
using ..VerticalGrids: stretched_vertical_cell_interfaces

#####
##### Bathymetry
#####

struct NeverworldBathymetry{G, B, C} <: Function
    grid :: G
    basin_depth_spline :: B
    channel_depth_spline :: C
    rim_width :: Float64
    slope_width :: Float64
    shelf_width :: Float64
    shelf_depth :: Float64
    abyssal_depth :: Float64
    southern_channel_boundary :: Float64
    northern_channel_boundary :: Float64
end

function NeverworldBathymetry(grid;
                              abyssal_depth = 4000,
                              shelf_depth = 200,
                              shelf_width = 2.5,
                              rim_width = shelf_width / 8,
                              slope_width = shelf_width,
                              southern_channel_boundary = -60,
                              northern_channel_boundary = -40)

    # Use grid spacing for "beach width"
    Δ = max(grid.Δλᶠᵃᵃ, grid.Δφᵃᶠᵃ)

    # Construct a cubic spline of the form `basin_depth(r)`, representing
    # the "basin component" of the Neverworld bathymetry where `r` is the distance
    # to the edge of the Neverworld (with units of degrees).
    r_coast = Δ
    r_beach = Δ + rim_width
    r_mid_shelf = Δ + rim_width + shelf_width/2
    r_shelf = Δ + rim_width + shelf_width
    r_abyss = Δ + rim_width + shelf_width + slope_width

    Nx, Ny, Nz = size(grid)
    x_max = xnode(f, c, c, Nx+1, 1, 1, grid) - xnode(f, c, c, 1, 1, 1, grid)
    y_max = ynode(c, f, c, 1, Ny+1, 1, grid) - ynode(c, f, c, 1, 1, 1, grid)
    r_max = max(x_max, y_max)

    basin_rim_distances = [0, r_coast,     r_beach,  r_mid_shelf,     r_shelf,       r_abyss,         r_max]
    basin_depths        = [0, 0,       shelf_depth,  shelf_depth, shelf_depth, abyssal_depth, abyssal_depth]
    basin_depth_spline = CubicSplineFunction{:x}(basin_rim_distances, basin_depths)

    # Construct cubic spline for the channel component using the "channel coordinate" δ.
    # Compared to the basin geometry, we omit the factor `Δ` representing the basin coastline.
    δ_beach = rim_width
    δ_shelf = rim_width + shelf_width
    δ_abyss = rim_width + shelf_width + slope_width
    δ_max = northern_channel_boundary - southern_channel_boundary
    channel_edge_distances = [0, δ_beach,       δ_shelf,       δ_abyss,         δ_max]
    channel_depths         = [0, rim_width, shelf_depth, abyssal_depth, abyssal_depth]
    channel_depth_spline = CubicSplineFunction{:y}(channel_edge_distances, channel_depths)

    return NeverworldBathymetry(grid,
                                basin_depth_spline,
                                channel_depth_spline,
                                Float64(rim_width),
                                Float64(slope_width),
                                Float64(shelf_width),
                                Float64(shelf_depth),
                                Float64(abyssal_depth),
                                Float64(southern_channel_boundary),
                                Float64(northern_channel_boundary))
end

const c = Center()
const f = Face()

function (nb::NeverworldBathymetry)(λ, φ)
    grid = nb.grid
    Nx, Ny, Nz = size(grid)

    # Four corners of the Neverworld
    λe = xnode(f, c, c, 1,       1, 1, grid)
    λw = xnode(f, c, c, Nx+1,    1, 1, grid)
    φs = ynode(c, f, c, 1,       1, 1, grid)
    φn = ynode(c, f, c, 1,    Ny+1, 1, grid)

    # Distance to the rim of the Neverworld
    r = min(λ - λe, λw - λ, φ - φs, φn - φ)
    basin_bottom_height = - nb.basin_depth_spline(r)

    # Channel coordinate: +0 inside, -0 outside
    δ = min(φ - nb.southern_channel_boundary, nb.northern_channel_boundary - φ)
    δ = max(0, δ)
    channel_bottom_height = - nb.channel_depth_spline(δ)

    # Intersect basin and channel
    bottom_height = min(basin_bottom_height, channel_bottom_height)

    return bottom_height
end

#####
##### Default vertical grid
#####

default_z = stretched_vertical_cell_interfaces(surface_layer_Δz = 5,
                                               surface_layer_height = 100,
                                               stretching_exponent = 1.02,
                                               minimum_depth = 4000)

#####
##### Boundary conditions
#####

# Default Neverworld wind stress profile
latitudes      = [-70,   -45,   -15,      0,    15,    45,  70]
zonal_stresses = [ +0,  -0.2,  +0.1,  +0.02,  +0.1,  -0.1,  +0] .* 1e-3 # kinematic wind stress
default_zonal_wind_stress = CubicSplineFunction{:y}(latitudes, zonal_stresses)

@inline cosine_target_buoyancy_distribution(φ, p) = p.Δb * cos(π * φ / p.Δφ)

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
                               #horizontal_resolution = 1, # degree
                               horizontal_size = (60, 70),
                               latitude = (-70, 0),
                               longitude = (0, 60),
                               z = default_z,
                               grid = nothing,
                               momentum_advection = VectorInvariant(vorticity_scheme   = WENO(),
                                                                    divergence_scheme  = WENO(),
                                                                    vertical_scheme    = WENO()),
                               tracer_advection = WENO(),
                               closure = CATKEVerticalDiffusivity(),
                               tracers = (:b, :e),
                               buoyancy = BuoyancyTracer(),
                               buoyancy_relaxation_time_scale = 4days,
                               target_buoyancy_distribution = cosine_target_buoyancy_distribution,
                               bottom_drag_coefficient = 2e-3,
                               equator_pole_buoyancy_difference = 0.06,
                               surface_buoyancy_gradient = 1e-5,
                               stratification_scale_height = 1000, # meters
                               time_step = 5minutes,
                               bathymetry = nothing,
                               zonal_wind_stress = default_zonal_wind_stress)


    if isnothing(grid)
        Nλ, Nφ = horizontal_size
        #Nλ = ceil(Int, longitude[2] - longitude[1] / horizontal_resolution)
        #Nφ = ceil(Int, latitude[2] - latitude[1] / horizontal_resolution)
        size = (Nλ, Nφ, length(z)-1)
        underlying_grid = LatitudeLongitudeGrid(arch; size, latitude, longitude, z, halo=(5, 5, 5))

        bathymetry = NeverworldBathymetry(underlying_grid)
        grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))
    end

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

