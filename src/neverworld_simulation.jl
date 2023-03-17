# Neverworld
#
# Ingredients:
#
#   * Zonally-periodic domain with continental shelves on all boundaries except Southern Ocean
#       * longitude = (0, 60)
#       * 4 configurations in latitude
#           - No Weddell Sea, half basin: latitude = (-60, 0)
#           - No Weddell Sea, full basin: latitude = (-60, 60)
#           - With Weddell Sea, half basin: latitude = (-70, 0)
#           - With Weddell Sea, full basin: latitude = (-70, 70)
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

using Oceananigans.Utils
using Oceananigans.Grids: node, halo_size
using Oceananigans.TurbulenceClosures: FluxTapering
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Operators: Δx, Δy, Az 
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization, ExplicitTimeDiscretization
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme

using Oceananigans.Fields: interpolate
using Oceananigans.Grids: xnode, ynode, halo_size

"""
    function z_faces_exp(; Nz = 69, Lz = 4000.0, e_folding = 0.06704463421863584)

generates an array of exponential z faces 

"""
function z_faces_exp(; Nz = 69, Lz = 4000.0, e_folding = 0.06704463421863584)
    z_faces   = zeros(Nz + 1)
    Nconstant = 11

    z_faces[1:Nconstant] .= 0:5:50

    for i in 1:(Nz + 1 - Nconstant)
        z_faces[i + Nconstant] = z_faces[i - 1 + Nconstant] + 5 * exp(e_folding * i)
    end

    z_faces    = - reverse(z_faces)
    z_faces[1] = - Lz

    return z_faces
end

"""
    function NeverworldGrid(arch, degree, FT::DataType = Float64; H = 5, longitude = (-2, 62), latitude = (-70, 0), bathymetry = bathymetry_without_ridge, longitudinal_extent = 60) 

builds a `LatitudeLongitudeGrid` with a specified `bathymetry`

Arguments
=========

- `arch` : architecture of the grid, can be `CPU()` or `GPU()`
- `resolution` : resolution in degrees.
- `FT` : (optional) floating point precision (default = `Float64`)

Keyword Arguments
=================

- `H` : halo size, `Int`
- `longitudinal_extent` : size of the actual domain in longitudinal direction, `Number`
- `longitude` : longitudinal extremes of the domain, `Tuple`. Note: this keyword must be at least `longitude_extent + resolution * 2H`
                to allow for correct advection stencils 
- `latitude` : latitudinal extremes of the domain
- `bathymetry` : function of `(λ, φ)` specifying the bottom height. Two bathymetry functions are implemented already:
                 `bathymetry_without_ridge` and `bathymetry_with_ridge`
- `z_faces` : array containing the z faces

"""
function NeverworldGrid(arch, resolution, FT::DataType = Float64; 
                        H = 5, longitudinal_extent = 60, 
                        longitude = (-2, 62), 
                        latitude = (-70, 0), 
                        bathymetry = bathymetry_without_ridge,
                        z_faces = z_faces_exp()) 

    Nx = Int((longitude[2] - longitude[1]) / resolution)
    Ny = Int((latitude[2]  - latitude[1]) / resolution)
    Nz = length(z_faces) - 1

    underlying_grid = LatitudeLongitudeGrid(arch, FT; size = (Nx, Ny, Nz),
                                            latitude,
                                            longitude,
                                            halo = (H, H, H),
                                            topology = (Periodic, Bounded, Bounded),
                                            z = z_faces)

    λ_grid = underlying_grid.λᶜᵃᵃ[1:Nx]
    φ_grid = underlying_grid.φᵃᶜᵃ[1:Ny]

    bathy = zeros(Nx, Ny)
    for (i, λ) in enumerate(λ_grid), (j, φ) in enumerate(φ_grid)
        bathy[i, j] = bathymetry(λ, φ; longitudinal_extent)
    end

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathy))
end

const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 
const ΔT   = 30.0
const fact = 5.0

"""
    function zonal_wind_stress(y, mid_wind)

returns the zonal wind as per https://egusphere.copernicus.org/preprints/2022/egusphere-2022-186/egusphere-2022-186.pdf
as a function of latitude `y` with  `mid_wind` the wind at the equator (`y = 0.0`)
    
"""
@inline function zonal_wind_stress(y, mid_wind)
    if y < -45
        return cubic_profile(y, -70.0, -45.0, 0.0, 0.2, 0.0, 0.0)
    elseif y < -15
        return cubic_profile(y, -45.0, -15.0, 0.2, -0.1, 0.0, 0.0)
    elseif y < 0
        return cubic_profile(y, -15.0, 0.0, -0.1, mid_wind, 0.0, 0.0)
    elseif y < 15
        return cubic_profile(y, 0.0, 15.0, mid_wind, -0.1, 0.0, 0.0)
    elseif y < 45
        return cubic_profile(y, 15.0, 45.0, -0.1, 0.1, 0.0, 0.0)
    else
        return cubic_profile(y, 45.0, 70.0, 0.1, 0.0, 0.0, 0.0)
    end
end

"""
    function salinity_flux(y, mid_flux)

returns the salinity flux as a function of latitude `y` 
(similar to https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020gl089135)
"""
@inline function salinity_flux(y)
    if y < -20
        return cubic_profile(y, -70.0, -20.0, -2e-8, 2e-8, 0.0, 0.0) .* 35.0
    elseif y < 0
        return cubic_profile(y, -20.0, 0.0, 2e-8, -4e-8, 0.0, 0.0) .* 35.0
    elseif y < 20
        return cubic_profile(y, 0.0, 20.0, -4e-8, 2e-8, 0.0, 0.0) .* 35.0
    else
        return cubic_profile(y, 20.0, 70.0, 2e-8, -2e-8, 0.0, 0.0) .* 35.0
    end
end

#####
##### Functions specifying initial conditions
#####

@inline exponential_profile(z; Δ = ΔB, Lz = Lz, h = h) = ( Δ * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) )

@inline parabolic_scaling(y) = - 1 / 70^2 * y^2 + 1
@inline atan_scaling(y)      = (atan(fact*((Ly + y)/Ly - 0.5)) / atan(fact * 0.5) + 1) /2

@inline initial_buoyancy_tangent(x, y, z)  = exponential_profile(z) * atan_scaling(y)
@inline initial_buoyancy_parabola(x, y, z) = exponential_profile(z) * parabolic_scaling(y) 

@inline initial_temperature_parabola(x, y, z) = exponential_profile(z; Δ = ΔT) * parabolic_scaling(y)

@inline function initial_salinity(y, mid_salinity)
    if y < -20
        return cubic_profile(y, -70.0, -20.0, 34.0, 37.0, 0.0, 0.0)
    elseif y < 0
        return cubic_profile(y, -20.0, 0.0, 37.0, 35.0, 0.0, 0.0)
    elseif y < 20
        return cubic_profile(y, 0.0, 20.0, 35.0, 37.0, 0.0, 0.0)
    else
        return cubic_profile(y, 20.0, 70.0, 37.0, 34.0, 0.0, 0.0)
    end
end

#####
##### Bottom drag boundary conditions
#####

@inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u))^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

#####
##### Wind stress boundary condition
#####

@inline surface_wind_stress(i, j, grid, clock, fields, τ) = τ[j]

#####
##### Tracers boundary condition
#####

@inline function buoyancy_top_relaxation(i, j, grid, clock, fields, p) 

    b = fields.b[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)

    return @inbounds p.λ * (b - p.initial_buoyancy(x, y, z))
end

@inline function temperature_top_relaxation(i, j, grid, clock, fields, p) 

    T  = fields.T[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)
    Trestoring = p.initial_temperature(x, y, z)

    return @inbounds p.λ * (T - Trestoring)
end

@inline function salinity_top_relaxation(i, j, grid, clock, fields, p) 

    S  = fields.S[i, j, grid.Nz]
    Srestoring = p.Ss[j]
    Sflux      = p.Fs[j]

    return @inbounds p.λ * (S - Srestoring) - Sflux
end

#####
##### Utility to concretize a boundary function into an array
#####

@inline function grid_specific_array(boundary_function, grid; scaling = 1000.0)

    Ny   = size(grid, 2)
    arch = architecture(grid)
    
    φ_grid = grid.φᵃᶜᵃ[1:Ny]

    τw = zeros(Ny)
    for (j, φ) in enumerate(φ_grid)
        τw[j] = boundary_function(φ, 0.0) ./ scaling
    end

    return arch_array(arch, -τw)
end

# The bathymetry is defined for a latitude range of -70 ≤ φ ≤ 0
# and a longitude range of 0 ≤ λ ≤ 60

""" 
    function cubic_profile(x, x1, x2, y1, y2, d1, d2)

returns a cubic function between points `(x1, y1)` and `(x2, y2)`
with derivative `d1` and `d2`
"""
@inline function cubic_profile(x, x1, x2, y1, y2, d1, d2)
    A = [x1^3  x1^2 x1 1
         x2^3  x2^2 x2 1
         3x1^2 2x1   1 0
         3x2^2 2x2   1 0]
          
    b = [y1, y2, d1, d2]

    coeff = A \ b

    return coeff[1] * x^3 + coeff[2] * x^2 + coeff[3] * x + coeff[4]
end

""" coastal rigde in longitude """
function coastal_ridge_x(x) 
    if x < 0.5
        return 0.0
    elseif x < 2.5 
        return -200.0
    elseif x < 5.0
        return cubic_profile(x, 2.5, 5.0, -200.0, -4000, 0.0, 0.0)
    else        
        return -4000.0 
    end
end

""" a sharp coast without a ridge """
function sharp_coast_x(x) 
    if x < 0.5
        return 0.0
    else        
        return -4000.0 
    end
end

""" coastal rigde in latitude """
function coastal_ridge_y(x) 
    if x < 0.5
        return cubic_profile(x, 0.0, 0.5, 0.0, -200, 0.0, 0.0)
    elseif x < 2.5 
        return -200.0
    elseif x < 5.0
        return cubic_profile(x, 2.5, 5.0, -200.0, -4000, 0.0, 0.0)
    else        
        return -4000.0 
    end
end

""" bottom atlantic rigde """
function bottom_ridge_x(x)
    if x < 20
        return -4000
    elseif x < 29
        return cubic_profile(x, 20.0, 29.0, -4000, -2000, 0.0, 0.0)
    elseif x < 31
        return -2000.0
    else 
        return -4000.0
    end
end

""" smoothed coasts for the inlet and outlet of the channel """
function bottom_ridge_xy(x, y)
    if y > - 30
        return bottom_ridge_x(x)
    elseif y > -50
        return cubic_profile(y, -30, -50, bottom_ridge_x(x), -4000, 0.0, 0.0)
    else
        return -4000.0
    end
end
        
""" scotia arc ridge """
function scotia_arc(x, y)
    radius = sqrt(x^2 + (y + 50)^2)
    if radius < 8
        return -4000.0
    elseif radius < 9
        return cubic_profile(radius, 8.0, 9.0, -4000.0, -2000.0, 0.0, 0.0)
    elseif radius < 11
        return -2000.0
    elseif radius < 12
        return cubic_profile(radius, 11.0, 12.0, -2000.0, -4000.0, 0.0, 0.0)
    else
        return -4000.0
    end
end

""" bathymetry without the atlantic ridge """
function bathymetry_without_ridge(x, y; longitudinal_extent = 60) 
    if x < 5 || x > 55
        if x < 0 
           x = 0.0
        end
        if x > 60
           x = 60.0
        end
        if y > -59 && y < -41 
            return  max(scotia_arc(x, y), 
                       coastal_ridge_x(sqrt(x^2 + (y + 59)^2)),
                       coastal_ridge_x(sqrt(x^2 + (y + 41)^2)), 
                       coastal_ridge_x(sqrt((longitudinal_extent - x)^2 + (y + 59)^2)),
                       coastal_ridge_x(sqrt((longitudinal_extent - x)^2 + (y + 41)^2)))
        else
            return max(coastal_ridge_x(x), 
                       coastal_ridge_y(70 + y),
                       coastal_ridge_y(70 + y),
                       coastal_ridge_x(longitudinal_extent - x), 
                       scotia_arc(x, y))
        end
    else
        return max(coastal_ridge_x(x),  
                   coastal_ridge_y(70 + y),
                   coastal_ridge_x(longitudinal_extent - x), 
                   scotia_arc(x, y))
    end
end

""" bathymetry with the atlantic ridge """
function bathymetry_with_ridge(x, y; longitudinal_extent = 60) 
    if x < 5 || x > 55
        if x < 0 
           x = 0.0
        end
        if x > 60
           x = 60.0
        end
        if y > -59 && y < -41 
            return  max(scotia_arc(x, y), 
                       coastal_ridge_x(sqrt(x^2 + (y + 59)^2)),
                       coastal_ridge_x(sqrt(x^2 + (y + 41)^2)), 
                       coastal_ridge_x(sqrt((longitudinal_extent - x)^2 + (y + 59)^2)),
                       coastal_ridge_x(sqrt((longitudinal_extent - x)^2 + (y + 41)^2)))
        else
            return max(coastal_ridge_x(x), 
                       coastal_ridge_x(longitudinal_extent - x), 
                       coastal_ridge_y(70 + y),
                       coastal_ridge_y(70 - y),
                       bottom_ridge_xy(x, y), 
                       bottom_ridge_xy(longitudinal_extent - x, y), 
                       scotia_arc(x, y))
        end
    else
        return max(coastal_ridge_x(x), 
                   coastal_ridge_x(longitudinal_extent - x), 
                   coastal_ridge_y(70 + y),
                   coastal_ridge_y(70 - y),
                   bottom_ridge_xy(x, y), 
                   bottom_ridge_xy(longitudinal_extent - x, y), 
                   scotia_arc(x, y))
    end
end

#####
##### Default parameterizations for the Neverworld simulation
#####

default_convective_adjustment  = RiBasedVerticalDiffusivity()
seawater_convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 0.2)
default_biharmonic_viscosity   = HorizontalDivergenceScalarBiharmonicDiffusivity(ν = geometric_νhb, discrete_form = true, parameters = 5days)
default_vertical_diffusivity   = VerticalScalarDiffusivity(ExplicitTimeDiscretization(), ν=1e-4, κ=1e-5)
default_slope_limiter          = FluxTapering(1e-2)

#####
##### Default momentum advection
#####

@inline upwind_vector_invariant(grid) = VectorInvariant(vorticity_scheme = WENO(), 
                                                       divergence_scheme = WENO(),
                                                         vertical_scheme = WENO(grid.underlying_grid))


"""
    function initialize_model!(model, Val(interpolate), initial_buoyancy, grid, orig_grid, init_file, buoyancymodel)

initializes the model according to 
1. interpolate or not on a finer/coarser grid `Val(interpolate)`
2. either `b` or `T` and `S`

"""
@inline initialize_model!(model, ::Val{false}, initial_buoyancy, grid, orig_grid, init_file, ::BuoyancyTracer) = set!(model, b = initial_buoyancy)

@inline function initialize_model!(model, ::Val{true}, initial_buoyancy, grid, orig_grid, init_file, ::BuoyancyTracer)
    Hx, Hy, Hz = halo_size(orig_grid)

    b_init = jldopen(init_file)["b/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    u_init = jldopen(init_file)["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    v_init = jldopen(init_file)["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    w_init = jldopen(init_file)["w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    
    @info "interpolating fields"
    b_init = interpolate_per_level(b_init, orig_grid, grid, (Center, Center, Center))
    u_init = interpolate_per_level(u_init, orig_grid, grid, (Face, Center, Center))
    v_init = interpolate_per_level(v_init, orig_grid, grid, (Center, Face, Center))
    w_init = interpolate_per_level(w_init, orig_grid, grid, (Center, Center, Face))

    set!(model, b=b_init, u=u_init, v=v_init, w=w_init) 
end

@inline initialize_model!(model, ::Val{false}, initial_profiles, grid, orig_grid, init_file, ::SeawaterBuoyancy) = set!(model, T = initial_profiles[1],  S = 35.0)

@inline function initialize_model!(model, ::Val{true}, initial_temperature, grid, orig_grid, init_file, ::SeawaterBuoyancy)
    Hx, Hy, Hz = halo_size(orig_grid)

    T_init = jldopen(init_file)["T/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    S_init = jldopen(init_file)["S/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    u_init = jldopen(init_file)["u/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    v_init = jldopen(init_file)["v/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    w_init = jldopen(init_file)["w/data"][Hx+1:end-Hx, Hy+1:end-Hy, Hz+1:end-Hz]
    
    @info "interpolating fields"
    T_init = interpolate_per_level(T_init, orig_grid, grid, (Center, Center, Center))
    S_init = interpolate_per_level(S_init, orig_grid, grid, (Center, Center, Center))
    u_init = interpolate_per_level(u_init, orig_grid, grid, (Face, Center, Center))
    v_init = interpolate_per_level(v_init, orig_grid, grid, (Center, Face, Center))
    w_init = interpolate_per_level(w_init, orig_grid, grid, (Center, Center, Face))

    set!(model, b=b_init, u=u_init, v=v_init, w=w_init) 
end

#=
# Calculate barotropic substeps based on barotropic CFL number and wave speed
function barotropic_substeps(Δt, grid, gravitational_acceleration; CFL = 0.7)
    wave_speed = sqrt(gravitational_acceleration * grid.Lz)
    local_Δ    = 1 / sqrt(1 / min_Δx(grid)^2 + 1 / min_Δy(grid)^2)

    return Int(ceil(2 * Δt / (CFL / wave_speed * local_Δ)))
end

free_surface = SplitExplicitFreeSurface(; substeps = barotropic_substeps(10minutes, grid, g_Earth))
=#

"""
    function weno_neverworld_simulation(; grid, 
                                          orig_grid = grid,
                                          drag_coefficient = 0.001,  
                                          buoyancy_relaxation_time_scale = 7days,
                                          convective_adjustment = default_convective_adjustment,
                                          biharmonic_viscosity  = default_biharmonic_viscosity,
                                          vertical_diffusivity  = default_vertical_diffusivity,
                                          gm_redi_diffusivities = nothing,
                                          tapering = default_slope_limiter,
                                          coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConservingScheme()),
                                          free_surface = ImplicitFreeSurface(),
                                          momentum_advection = upwind_vector_invariant(grid),
    				                      tracer_advection   = WENO(grid), 
                                          interp_init = false,
                                          init_file = nothing,
                                          time_step = 5minutes,
                                          stop_time = 10years,
                                          buoyancy_boundary_conditions = true,
                                          velocity_boundary_conditions = true,
                                          initial_buoyancy = initial_buoyancy_tangent,
    				                      wind_stress  = zonal_wind_stress,
                                          tracers = :b
                                          )

returns a neverworld `simulation` using a `BuoyancyTracer` buoyancy model

Keyword arguments
=================

- `grid` : the neverworld `LatitudeLongitudeGrid`
- `orig_grid` : if `interp_init == true` we need the original grid to perform the interpolation
- `μ_drag` : drag coefficient (default 1e-3)
- `λ_buoy` : restoring time for buoyancy restoration (default = `7days`)
- `convective_adjustment` : boundary layer parameterization (default = `RiBasedVerticalDiffusivity()`)
- `biharmonic_viscosity` : horizontal momentum diffusion
- `vertical_diffusivity` : vertical (background) momentum and tracer diffusion 
- `gm_redi_diffusivities` : a tuple containing `(κᴳᴹ, κᴿᴱᴰᴵ)`, if `nothing` no gm closure is used 
                            (default = `nothing`)
- `tapering` : tapering for the gm-redi closure (default = `FluxTapering(1e-2)`)
- `coriolis` : coriolis scheme
- `free_surface` : free surface scheme
- `momentum_advection` : momentum advection scheme
- `tracer_advection` : tracer advection scheme
- `interp_init` : if `false` the simulation will be initialized with `initial_buoyancy` and quiescent flow,
                  if `true` the simulation will be initialized by interpolating the values in `init_file`
- `init_file` : initial file used to interpolate variables
- `Δt` : time step size (constant)
- `stop_time` : stop time of the simulation
- `buoyancy_boundary_conditions` : if `true` we use restoring boundary conditions for buoyancy at the top,
                                   otherwise we use zero flux
- `velocity_boundary_conditions` : if `true` we impose a wind stress and drag boundary conditions,
                                   otherwise we use zero flux
- `initial_buoyancy` : initial buoyancy profile. The top layer is used as a reference for the restoring boundary condition
- `wind_stress` : wind stress profile. It is used as a top boundary condition for `u`
- `tracers` : tracers in the simulation
"""
function weno_neverworld_simulation(; grid, 
                                      orig_grid = grid,
                                      μ_drag = 0.001,  
                                      λ_buoy = 7days,
                                      convective_adjustment = default_convective_adjustment,
                                      biharmonic_viscosity  = default_biharmonic_viscosity,
                                      vertical_diffusivity  = default_vertical_diffusivity,
                                      gm_redi_diffusivities = nothing,
                                      tapering = default_slope_limiter,
                                      coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConservingScheme()),
                                      free_surface = ImplicitFreeSurface(),
                                      momentum_advection = upwind_vector_invariant(grid),
				                      tracer_advection   = WENO(grid.underlying_grid), 
                                      interp_init = false,
                                      init_file = nothing,
                                      Δt = 5minutes,
                                      stop_time = 10years,
                                      buoyancy_boundary_conditions = true,
                                      velocity_boundary_conditions = true,
                                      initial_buoyancy = initial_buoyancy_tangent,
				                      wind_stress  = zonal_wind_stress,
                                      tracers = :b
                                      )

    # Initializing boundary conditions

    @info "specifying boundary conditions..."

    @apply_regionally τw = grid_specific_array(wind_stress, grid)

    u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τw)

    # Quadratic bottom drag:

    if velocity_boundary_conditions
        drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)
        drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = μ_drag)

        u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u) 
        v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v) 

        u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ_drag)
        v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ_drag)

        u_bcs = FieldBoundaryConditions(bottom = u_bottom_drag_bc, immersed = u_immersed_bc, top = u_wind_stress_bc)
        v_bcs = FieldBoundaryConditions(bottom = v_bottom_drag_bc, immersed = v_immersed_bc)
    else
        u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
        v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    end

    Δz_top = CUDA.@allowscalar Δzᶜᶜᶜ(1, 1, grid.Nz, grid)
    v_pump = Δz_top / λ_buoy

    b_top_relaxation_bc = FluxBoundaryCondition(buoyancy_top_relaxation, discrete_form=true, parameters = (; λ = v_pump, initial_buoyancy))

    if buoyancy_boundary_conditions
        b_bcs = FieldBoundaryConditions(top = b_top_relaxation_bc)
    else
        b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0))
    end

    #####
    ##### Closures
    #####

    @info "specifying closures..."

    if gm_redi_diffusivities isa Nothing
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment)
    else
        κᴳ, κᴿ = gm_redi_diffusivities
        isopycnal_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew = κᴳ, κ_symmetric = κᴿ, slope_limiter = tapering)
        closure = (vertical_diffusivity, biharmonic_viscosity, convective_adjustment, isopycnal_closure)
    end

    #####
    ##### Model setup
    #####

    @info "building model..."            

    model = HydrostaticFreeSurfaceModel(; grid, free_surface, coriolis, closure, tracers, momentum_advection, tracer_advection, 
                                          boundary_conditions = (; u = u_bcs, v = v_bcs, b = b_bcs), 
                                          buoyancy = BuoyancyTracer())

    #####
    ##### Model initialization
    #####

    @info "initializing prognostic variables from $(interp_init ? init_file : "scratch")"
    initialize_model!(model, Val(interp_init), initial_buoyancy, grid, orig_grid, init_file, BuoyancyTracer())

    simulation = Simulation(model; Δt, stop_time)

    @show start_time = [time_ns()]

    function progress(sim)
        sim.model.clock.iteration == 1

        wall_time = (time_ns() - start_time[1]) * 1e-9

        u, v, w = sim.model.velocities

        @info @sprintf("Time: % 12s, it: %d, max(|u|, |v|, |w|): (%.2e, %.2e , %.2e) ms⁻¹, Δt: %.2e s, wall time: %s", 
            prettytime(sim.model.clock.time),
	    sim.model.clock.iteration, maximum(abs, u), maximum(abs, v), maximum(abs, w), sim.Δt,
            prettytime(wall_time))

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

    return simulation
end

"""
    function run_simulation!(simulation; interp_init = false, init_file = nothing) 

runs the `simulation`. if `interp_init` is `false` and `init_file` is not a `Nothing`
we pickup from `init_file`
"""
function run_simulation!(simulation; interp_init = false, init_file = nothing) 
    
    init = interp_init ? true : (init_file isa Nothing ? true : false)

    Δt    = simulation.Δt
    model = simulation.model 
        
    if init
        @info "running simulation from zero-velocity initial conditions"
        run!(simulation)
    else
        @info "running simulation from $init_file"
        update_simulation_clock!(simulation, init_file)
        run!(simulation, pickup=init_file)
    end
    
    @info """
        Simulation took $(prettytime(simulation.run_wall_time))
        Free surface: $(typeof(model.free_surface).name.wrapper)
        Time step: $(prettytime(Δt))
    """
end
