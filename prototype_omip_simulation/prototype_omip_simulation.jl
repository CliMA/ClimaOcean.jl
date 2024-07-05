using Printf
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Operators
using Oceananigans.Operators: ℑxyz
using Oceananigans: architecture
using Oceananigans.Grids: on_architecture, znode
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving, fᶠᶠᵃ
using Oceananigans.BuoyancyModels: ∂y_b, ∂z_b
using Oceananigans.Units
using ClimaOcean
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: Radiation, SimilarityTheoryTurbulentFluxes
using ClimaOcean.VerticalGrids: exponential_z_faces
using ClimaOcean.JRA55
using ClimaOcean.ECCO
using ClimaOcean.JRA55: JRA55NetCDFBackend, JRA55_prescribed_atmosphere
using ClimaOcean.ECCO: ECCO_restoring_forcing, ECCO4Monthly, ECCO2Daily, ECCOMetadata
using ClimaOcean.Bathymetry
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: LatitudeDependentAlbedo
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

import ClimaOcean: stateindex

using CFTime
using Dates

include("tripolar_specific_methods.jl")

#####
##### Global Ocean at 1/6th of a degree
#####

bathymetry_file = nothing # "bathymetry_tmp.jld2"

# 60 vertical levels
z_faces = exponential_z_faces(Nz=30, depth=6000)

Nx = 360
Ny = 180
Nz = length(z_faces) - 1

arch = CPU() #Distributed(GPU(), partition = Partition(2))

grid = TripolarGrid(arch; 
                    size = (Nx, Ny, Nz), 
                    halo = (7, 7, 7), 
                    z = z_faces, 
                    north_poles_latitude = 55,
                    first_pole_longitude = 75)

bottom_height = retrieve_bathymetry(grid, bathymetry_file; 
                                    minimum_depth = 10,
                                    dir = "./",
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 0)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 

#####
##### The Ocean component
#####                             

free_surface = SplitExplicitFreeSurface(grid; substeps = 75)

#####
##### Add restoring to ECCO fields for temperature and salinity in the artic and antarctic
#####

# Build a mask that goes from 0 to 1 as a cubic function of φ between
# 70 degrees and 90 degrees and zero derivatives at 70 and 90.
x₁ = 70
x₂ = 90
y₁ = 0
y₂ = 1

A⁺ = [ x₁^3   x₁^2  x₁ 1
       x₂^3   x₂^2  x₂ 1
       3*x₁^2 2*x₁  1  0
       3*x₂^2 2*x₂  1  0]
           
b⁺ = [y₁, y₂, 0, 0]
c⁺ = A⁺ \ b⁺

# Coefficients for the cubic mask
const c₁⁺ = c⁺[1]
const c₂⁺ = c⁺[2]
const c₃⁺ = c⁺[3]
const c₄⁺ = c⁺[4]

const c₁⁻ = - c⁺[1]
const c₂⁻ = c⁺[2]
const c₃⁻ = - c⁺[3]
const c₄⁻ = c⁺[4]

@inline mask_f(λ, φ, z) = ifelse(φ >=  70, c₁⁺ * φ^3 + c₂⁺ * φ^2 + c₃⁺ * φ + c₄⁺,
                          ifelse(φ <= -70, c₁⁻ * φ^3 + c₂⁻ * φ^2 + c₃⁻ * φ + c₄⁻, zero(eltype(φ))))

mask = CenterField(grid)
set!(mask, mask_f)

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(1993, 12, 1)

temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
salinity    = ECCOMetadata(:salinity,    dates, ECCO4Monthly())

FT = ECCO_restoring_forcing(temperature; mask, grid, architecture = arch, timescale = 30days)
FS = ECCO_restoring_forcing(salinity;    mask, grid, architecture = arch, timescale = 30days)

forcing = (; T = FT, S = FS)

#####
##### Adding custom closures!!!!
#####

# Vertical Closure!
vertical_closure = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 0.1,
                                                           background_κz = 1e-5,
                                                           convective_νz = 0.1,
                                                           background_νz = 1e-5)

buoyancy = SeawaterBuoyancy(; gravitational_acceleration = Oceananigans.BuoyancyModels.g_Earth, 
                              equation_of_state= TEOS10EquationOfState(; reference_density = 1020))

# Parameters for the κskew function
gm_parameters = (; max_C = 20, 
                   min_C = 0.23,
                   K₀ᴳᴹ  = 1e3,
                   buoyancy,
                   coriolis = HydrostaticSphericalCoriolis()
                 )

# Custom GM coefficient, function of `i, j, k, grid, clock, fields, parameters`
@inline function κskew(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, p) 
     
     β = ∂yᶜᶜᶜ(i, j, k, grid, fᶠᶠᵃ, p.coriolis)

    ∂ʸb = ℑxyz(i, j, k, grid, (Center(), Face(), Center()), (ℓx, ℓy, ℓz), ∂y_b, p.buoyancy, fields)
    ∂ᶻb = ℑxyz(i, j, k, grid, (Center(), Center(), Face()), (ℓx, ℓy, ℓz), ∂z_b, p.buoyancy, fields)

    Sʸ = ∂ʸb / ∂ᶻb 
    Sʸ = ifelse(isnan(Sʸ), zero(grid), Sʸ)

    z  = znode(k, grid.underlying_grid, Center())
    C  = min(max(p.min_C, 1 - β / Sʸ * z), p.max_C)

    return p.K₀ᴳᴹ * C
end

gerdes_koberle_willebrand_tapering = FluxTapering(1e-1)
horizontal_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew = κskew,
                                                       κ_symmetric = 1000,
                                                       skew_discrete_form = true,
                                                       skew_loc = (nothing, nothing, nothing),
                                                       parameters = gm_parameters,
                                                       slope_limiter = gerdes_koberle_willebrand_tapering,
                                                       required_halo_size = 3)


#####
##### Building and initializing the ocean simulation!
#####

closure = (vertical_closure, horizontal_closure)
ocean   = ocean_simulation(grid; 
                           tracers = (:T, :S, :e),
                           free_surface, 
                           forcing, 
                           closure) 

#####
##### The atmosphere
#####

backend    = JRA55NetCDFBackend(4) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation(arch)

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

wall_time = [time_ns()]

function progress(sim) 
    u, v, w = sim.model.velocities  
    T, S, e = sim.model.tracers

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(interior(u)), maximum(interior(v)), maximum(interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

ocean.callbacks[:progress] = Callback(progress, IterationInterval(10))

fluxes = (u = model.velocities.u.boundary_conditions.top.condition,
          v = model.velocities.v.boundary_conditions.top.condition,
          T = model.tracers.T.boundary_conditions.top.condition,
          S = model.tracers.S.boundary_conditions.top.condition)

ocean.output_writers[:fluxes] = JLD2OutputWriter(model, fluxes,
                                                  schedule = TimeInterval(0.5days),
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32},
                                                  filename = "surface_fluxes")

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                  schedule = TimeInterval(0.5days),
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32},
                                                  filename = "surface",
                                                  indices = (:, :, grid.Nz))

ocean.output_writers[:snapshots] = JLD2OutputWriter(model, merge(model.tracers, model.velocities),
                                                    schedule = TimeInterval(10days),
                                                    overwrite_existing = true,
                                                    array_type = Array{Float32},
                                                    filename = "snapshots")

ocean.output_writers[:checkpoint] = Checkpointer(model, 
                                                 schedule = TimeInterval(60days),
                                                 overwrite_existing = true,
                                                 prefix = "checkpoint")

restart = nothing
model   = ocean.model
initial_date = dates[1]

if isnothing(restart) # Warm up!
     set!(model, 
          T = ECCOMetadata(:temperature, initial_date, ECCO4Monthly()),
          S = ECCOMetadata(:salinity,    initial_date, ECCO4Monthly()),
          e = 1e-6)

     ocean.Δt = 10
     ocean.stop_iteration = 1
     wizard = TimeStepWizard(; cfl = 0.1, max_Δt = 90, max_change = 1.1)
     ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

     stop_time = 30days

     coupled_simulation = Simulation(coupled_model; Δt=1, stop_time)

     run!(coupled_simulation)
else
     set!(model, restart)
end

wizard = TimeStepWizard(; cfl = 0.3, max_Δt = 900, max_change = 1.1)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's reset the maximum number of iterations
coupled_model.ocean.stop_time = 7200days
coupled_simulation.stop_time = 7200days
coupled_model.ocean.stop_iteration = Inf
coupled_simulation.stop_iteration = Inf

run!(coupled_simulation)
