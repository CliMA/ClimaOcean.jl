function default_latlon_closure()
    catke = default_ocean_closure()
    horizontal_viscosity = HorizontalScalarDiffusivity(ν=5e4)
    vertical_diffusivity = VerticalScalarDiffusivity(ν=1e-5, κ=henyey_diffusivity)
    return (catke, horizontal_viscosity, vertical_diffusivity)
end

"""
    latitude_longitude_ocean(arch = CPU(); zstar=true, kwargs...)

Construct an ocean `Simulation` on a one-degree `LatitudeLongitudeGrid` (360×150)
spanning 75°S–75°N with realistic bathymetry.
"""
function latitude_longitude_ocean(arch = CPU();
                                  zstar = true,
                                  momentum_advection = VectorInvariant(),
                                  tracer_advection = WENO(order=7),
                                  closure = default_latlon_closure(),
                                  kwargs...)

    z = vertical_coordinate(; zstar)

    grid = LatitudeLongitudeGrid(arch;
                                 size = (360, 150, 60),
                                 z,
                                 halo = (7, 7, 7),
                                 latitude = (-75, 75),
                                 longitude = (0, 360))

    bottom_height = regrid_bathymetry(grid;
                                      minimum_depth = 10,
                                      interpolation_passes = 5,
                                      major_basins = 3)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height);
                                active_cells_map = true)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            closure,
                            kwargs...)
end
