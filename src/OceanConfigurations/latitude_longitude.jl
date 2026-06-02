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
                                  Nz = 60,
                                  depth = 6000,
                                  momentum_advection = VectorInvariant(),
                                  tracer_advection = WENO(order=7),
                                  closure = default_latlon_closure(),
                                  halo = (7, 7, 7),
                                  minimum_depth = 10,
                                  interpolation_passes = 5,
                                  z = nothing,
                                  additional_surface_fluxes = nothing,
                                  kwargs...)

    if isnothing(z)
        z = vertical_coordinate(; Nz, depth, zstar)
    end

    grid = LatitudeLongitudeGrid(arch;
                                 size = (360, 150, Nz),
                                 z,
                                 halo,
                                 latitude = (-75, 75),
                                 longitude = (0, 360))

    bottom_height = regrid_bathymetry(grid;
                                      minimum_depth,
                                      interpolation_passes,
                                      major_basins = 3)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height);
                                active_cells_map = true)

    asf = resolve_surface_fluxes(additional_surface_fluxes, arch, grid)
    flux_kw = isnothing(asf) ? (;) : (; additional_surface_fluxes = asf)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            closure,
                            flux_kw...,
                            kwargs...)
end
