using ClimaOcean
using Oceananigans

Nx, Ny, Nz = 32, 32, 10
z = exponential_z_faces(Nz=Nz, depth=6000, h=34)
grid = Oceananigans.OrthogonalSphericalShellGrids.TripolarGrid(CPU(); size=(Nx, Ny, Nz), halo=(7, 7, 7), z)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom((x, y) -> -5000))
ocean = ocean_simulation(grid)

@time @eval OceanSeaIceModel(ocean)

