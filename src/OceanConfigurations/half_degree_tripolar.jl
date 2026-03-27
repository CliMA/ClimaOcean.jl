function default_half_degree_closure()
    catke = default_ocean_closure()
    eddy  = IsopycnalSkewSymmetricDiffusivity(κ_skew=500, κ_symmetric=200)
    horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true, parameters=40days)
    vertical_diffusivity = VerticalScalarDiffusivity(ν=1e-5, κ=henyey_diffusivity)
    return (catke, eddy, horizontal_viscosity, vertical_diffusivity)
end

"""
    half_degree_tripolar_ocean(arch = CPU(); zstar=false, kwargs...)

Construct an ocean `Simulation` on a half-degree (720×360) `TripolarGrid` with
realistic bathymetry and production-ready closures (CATKE + Gent-McWilliams +
biharmonic viscosity).
"""
function half_degree_tripolar_ocean(arch = CPU();
                                    zstar = false,
                                    momentum_advection = WENOVectorInvariant(order=5),
                                    tracer_advection = WENO(order=7),
                                    closure = default_half_degree_closure(),
                                    kwargs...)

    z = vertical_coordinate(; zstar)

    grid = TripolarGrid(arch;
                        size = (720, 360, 60),
                        z,
                        halo = (7, 7, 7))

    bottom_height = regrid_bathymetry(grid;
                                      minimum_depth = 20,
                                      major_basins = 1,
                                      interpolation_passes = 25)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height);
                                active_cells_map = true)

    free_surface = SplitExplicitFreeSurface(grid; substeps=150)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            free_surface,
                            closure,
                            kwargs...)
end
