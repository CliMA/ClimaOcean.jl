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

using ClimaOcean: CubicSplineFunction,
                  u_bottom_drag, u_immersed_bottom_drag,
                  v_bottom_drag, v_immersed_bottom_drag
                  
using ..VerticalGrids: stretched_vertical_cell_interfaces

const c = Center()
const f = Face()

#####
##### Geometry
#####

struct Point{T}
    x :: T
    y :: T
end

struct LineSegment{T}
    p₁ :: Point{T}
    p₂ :: Point{T}
end

struct Line{T}
    p₁ :: Point{T}
    p₂ :: Point{T}
end

distance(p₁::Point, p₂::Point) = sqrt((p₁.x - p₂.x)^2 + (p₁.y - p₂.y)^2)

"""
   distance(point::Point, linesegment::LineSegment)

Return the distance between a `point` and a `linesegment`, that is the shortest distance
 of the `point` with and of the points within the line segment.
"""
function distance(point::Point, linesegment::LineSegment)
    x, y = point.x, point.y
    x₁, y₁ = linesegment.p₁.x, linesegment.p₁.y
    x₂, y₂ = linesegment.p₂.x, linesegment.p₂.y

    # Line segment vector components
    ℓˣ = x₂ - x₁
    ℓʸ = y₂ - y₁

    # Fractional increment to closest segment point
    ϵ = ((x - x₁) * ℓˣ + (y - y₁) * ℓʸ) / (ℓˣ^2 + ℓʸ^2)
    ϵ = clamp(ϵ, 0, 1)

    # Closest segment point
    x′ = x₁ + ϵ * ℓˣ
    y′ = y₁ + ϵ * ℓʸ
    
    return distance(Point(x, y), Point(x′, y′))
end

function distance(point::Point, line::Line)
    x, y = point.x, point.y
    x₁, y₁ = line.p₁.x, line.p₁.y
    x₂, y₂ = line.p₂.x, line.p₂.y
    num = abs((x₂ - x₁) * (y₁ - y) - (y₂ - y₁) * (x₁ - x))
    den = sqrt((x₂ - x₁)^2 + (y₂ - y₁)^2)
    return num / den
end

#####
##### Bathymetry
#####

struct NeverworldBathymetry{G, B} <: Function
    grid :: G
    coastline_spline :: B
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
    r_mid_shelf = Δ + rim_width + shelf_width / 2
    r_shelf = Δ + rim_width + shelf_width
    r_abyss = Δ + rim_width + shelf_width + slope_width

    Nx, Ny, Nz = size(grid)
    x_max = xnode(Nx+1, 1, 1, grid, f, c, c) - xnode(1, 1, 1, grid, f, c, c)
    y_max = ynode(1, Ny+1, 1, grid, c, f, c) - ynode(1, 1, 1, grid, c, f, c)
    r_max = max(x_max, y_max)

    basin_rim_distances = [0, r_coast,     r_beach,  r_mid_shelf,     r_shelf,       r_abyss,         r_max]
    basin_depths        = [0, 0,       shelf_depth,  shelf_depth, shelf_depth, abyssal_depth, abyssal_depth]
    coastline_spline = CubicSplineFunction{:x}(basin_rim_distances, basin_depths)

    return NeverworldBathymetry(grid,
                                coastline_spline,
                                Float64(rim_width),
                                Float64(slope_width),
                                Float64(shelf_width),
                                Float64(shelf_depth),
                                Float64(abyssal_depth),
                                Float64(southern_channel_boundary),
                                Float64(northern_channel_boundary))
end

function (nb::NeverworldBathymetry)(λ, φ)
    grid = nb.grid
    Nx, Ny, Nz = size(grid)

    # Four corners of the Neverworld
    λw = xnode(1,       1, 1, grid, f, c, c)
    λe = xnode(Nx+1,    1, 1, grid, f, c, c)
    φs = ynode(1,       1, 1, grid, c, f, c)
    φn = ynode(1,    Ny+1, 1, grid, c, f, c)

    # Draw lines along the six coasts of the Neverworld
    northern_vertices = (Point(λw, nb.northern_channel_boundary),
                         Point(λw, φn),
                         Point(λe, φn),
                         Point(λe, nb.northern_channel_boundary))

    southern_vertices = (Point(λe, nb.southern_channel_boundary),
                         Point(λe, φs),
                         Point(λw, φs),
                         Point(λw, nb.southern_channel_boundary))

    coastlines = [LineSegment(northern_vertices[1], northern_vertices[2]),
                  LineSegment(northern_vertices[2], northern_vertices[3]),
                  LineSegment(northern_vertices[3], northern_vertices[4]),
                  LineSegment(southern_vertices[1], southern_vertices[2]),
                  LineSegment(southern_vertices[2], southern_vertices[3]),
                  LineSegment(southern_vertices[3], southern_vertices[4])]

    # Minimum distance to the six rims of the Neverworld
    p = Point(λ, φ)
    r = minimum(distance(p, coastline) for coastline in coastlines)
    
    bottom_height = - nb.coastline_spline(r)

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
    φ = ynode(i, j, k, grid, c, c, c)
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
        size = (Nλ, Nφ, length(z)-1)
        underlying_grid = LatitudeLongitudeGrid(arch; size, latitude, longitude, z, halo=(5, 5, 5))

        bathymetry = NeverworldBathymetry(underlying_grid)
        grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))
    end

    N² = surface_buoyancy_gradient
    h  = stratification_scale_height 
    Δb = equator_pole_buoyancy_difference 
    t★ = buoyancy_relaxation_time_scale 
    Δφ = abs(latitude[1])
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
                                        momentum_advection, tracer_advection, closure,
                                        boundary_conditions = (b=b_bcs, u=u_bcs, v=v_bcs))

    bᵢ(x, y, z) = Δb + N² * h * (exp(z / h) - 1)
    set!(model, b=bᵢ, e=1e-6)

    simulation = Simulation(model; Δt=time_step)

    return simulation
end

