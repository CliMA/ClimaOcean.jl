using Printf
using Oceananigans
using Oceananigans.Units
using ClimaOcean
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans: architecture
using ClimaOcean
using ClimaOcean.ECCO2
using Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.Coriolis: ActiveCellEnstrophyConserving
using Oceananigans.Units
using ClimaOcean.OceanSimulations
using ClimaOcean.OceanSeaIceModels
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: Radiation
using ClimaOcean.VerticalGrids: exponential_z_faces
using ClimaOcean.JRA55
using ClimaOcean.JRA55: JRA55NetCDFBackend, JRA55_prescribed_atmosphere
using ClimaOcean.Bathymetry

include("three_dimensional_interpolate_tripolar.jl")

#####
##### Global Ocean at 1/4th of a degree
#####

bathymetry_file = nothing # "bathymetry_tmp.jld2"

# 60 vertical levels
z_faces = exponential_z_faces(Nz=60, depth=6500)

Nx = 720
Ny = 300
Nz = length(z_faces) - 1

arch = CPU() #Distributed(GPU(), partition = Partition(2))

grid = TripolarGrid(arch; size = (Nx, Ny, Nz), halo = (7, 7, 7), z = z_faces)

bottom_height = retrieve_bathymetry(grid, bathymetry_file; 
                                    minimum_depth = 10,
                                    dir = "./",
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 0)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 

#####
##### The Ocean component
#####                             

const Lz = grid.Lz
const  h = Nz / 4.5

@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h))
@inline νz(x, y, z, t) = 1e-4 + (5e-3 - 1e-4) * exponential_profile(z; Lz, h)

free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt = 20minutes)
vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = 5e-5, ν = 5e-3)

closure = (RiBasedVerticalDiffusivity(), vertical_diffusivity)

ocean = ocean_simulation(grid; free_surface, closure) 
model = ocean.model

set!(model, 
     T = ECCO2Metadata(:temperature),
     S = ECCO2Metadata(:salinity))

#####
##### The atmosphere
#####

backend    = JRA55NetCDFBackend(10) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation  = Radiation()

sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()

coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

wall_time = [time_ns()]

function progress(sim) 
    u, v, w = sim.model.velocities  
    T, S = sim.model.tracers

    Tmax = maximum(abs, T)
    Tmin = minimum(T)
    umax = maximum(abs, u), maximum(abs, v), maximum(abs, w)
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end
