using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.OrthogonalSphericalShellGrids

r_faces = ClimaOcean.exponential_z_faces(; Nz=30, h=10, depth=3000)
z_faces = MutableVerticalDiscretization(r_faces)

Nx = 180 # longitudinal direction -> 250 points is about 1.5ᵒ resolution
Ny = 180 # meridional direction -> same thing, 48 points is about 1.5ᵒ resolution
Nz = length(r_faces) - 1

grid = RotatedLatitudeLongitudeGrid(size = (Nx, Ny, Nz), 
                                    latitude = (-45, 45),
                                    longitude = (-45, 45),
                                    z = z_faces,
                                    north_pole = (0, 0),
                                    topology = (Bounded, Bounded, Bounded))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

@inline x_domain(grid::RotatedLatitudeLongitudeGrid) = domain(topology(grid, 1)(), grid.Nx, grid.xᶠᵃᵃ)
@inline y_domain(grid::RotatedLatitudeLongitudeGrid) = domain(topology(grid, 2)(), grid.Ny, grid.yᵃᶠᵃ)

#####
##### Ocean model
#####

momentum_advection = WENOVectorInvariant(order=5) 
tracer_advection   = Centered()

free_surface = SplitExplicitFreeSurface(grid; substeps=30) 

ocean = ocean_simulation(grid; 
                         momentum_advection, 
                         tracer_advection, 
                         free_surface)

#####
##### Sea-ice model
#####

sea_ice = sea_ice_simulation(grid; dynamics=nothing, advection=nothing) 

#####
##### Atmosphere model
#####

atmosphere = JRA55PrescribedAtmosphere(; backend=JRA55NetCDFBackend(40))
