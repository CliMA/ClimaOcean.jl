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

include("correct_oceananigans.jl")

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
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 1)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map = true) 

#####
##### The Ocean component
#####                             

const Lz = grid.Lz
const  h = Nz / 4.5

@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h))
@inline νz(x, y, z, t) = 1e-4 + (5e-3 - 1e-4) * exponential_profile(z; Lz, h)

free_surface = SplitExplicitFreeSurface(; grid, cfl=0.7, fixed_Δt = 20minutes)
vertical_diffusivity = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = 5e-5, ν = 5e-3)

closure = (RiBasedVerticalDiffusivity(), vertical_diffusivity)

ocean = ocean_simulation(grid; free_surface, closure) 
model = ocean.model