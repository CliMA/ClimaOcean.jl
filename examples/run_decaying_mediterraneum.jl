using GLMakie
using Oceananigans
using ClimaOcean.Bathymetry: regrid_bathymetry
using ClimaOcean.ECCO2: ecco2_field, ecco2_center_mask
using ClimaOcean.VerticalGrids: stretched_vertical_faces, PowerLawStretching

using ClimaOcean.InitialConditions: three_dimensional_regrid!

#####
##### Regional Mediterranean grid 
#####

# A stretched vertical grid with a Δz of 1.5 meters in the first 50 meters
z = stretched_vertical_faces(minimum_depth = 5000, 
                             surface_layer_Δz = 1.75, 
                             stretching = PowerLawStretching(1.070), 
                             surface_layer_height = 50)

Nx = 20 * 55 # 1 / 20th of a degree
Ny = 20 * 25
Nz = length(z) - 1

grid = LatitudeLongitudeGrid(CPU();
                             size = (Nx, Ny, Nz),
                             latitude = (25, 50),
                             longitude = (-10, 45),
                             z,
                             halo = (7, 7, 7))

h = regrid_bathymetry(grid, height_above_water=1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(h))

# All this can be done in a `set_tracers_from_ecco!(model; )` 
# Download ecco tracer fields
Tecco = ecco2_field(:temperature);
Secco = ecco2_field(:salinity);

ecco_tracers = (; Tecco, Secco)

# Make sure all values are extended properly before regridding
adjust_tracers!(ecco_tracers; mask = ecco2_center_mask())

T = CenterField(grid);
S = CenterField(grid);

three_dimensional_regrid!(T, Tecco)
three_dimensional_regrid!(S, Secco)
