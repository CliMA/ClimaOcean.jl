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

include("correct_oceananigans.jl")

#####
##### Global Ocean at 1/12th of a degree
#####

bathymetry_file = nothing # "bathymetry_tmp.jld2"

# 100 vertical levels
z_faces = exponential_z_faces(Nz=60, depth=6500)

Nx = 720
Ny = 300
Nz = length(z_faces) - 1

arch = CPU() #Distributed(GPU(), partition = Partition(2))

underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), halo = (7, 7, 7))

bottom_height = retrieve_bathymetry(grid, bathymetry_file; 
                                    minimum_depth = 10,
                                    interpolation_passes = 20,
                                    connected_regions_allowed = 3)
 
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map = true) 