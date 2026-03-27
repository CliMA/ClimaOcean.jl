function default_one_degree_closure()
    catke = default_ocean_closure()
    eddy  = IsopycnalSkewSymmetricDiffusivity(κ_skew=500, κ_symmetric=200)
    horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true, parameters=15days)
    vertical_diffusivity = VerticalScalarDiffusivity(ν=1e-5, κ=henyey_diffusivity)
    return (catke, eddy, horizontal_viscosity, vertical_diffusivity)
end

"""
    one_degree_tripolar_ocean(arch = CPU(); zstar=false, kwargs...)

Construct an ocean `Simulation` on a one-degree (360×180) `TripolarGrid` with
realistic bathymetry and closures suitable for coarser resolution
(CATKE + Gent-McWilliams + biharmonic viscosity).
"""
function one_degree_tripolar_ocean(arch = CPU();
                                   zstar = false,
                                   momentum_advection = WENOVectorInvariant(order=5),
                                   tracer_advection = WENO(order=5),
                                   closure = default_one_degree_closure(),
                                   kwargs...)

    z = vertical_coordinate(; zstar)

    grid = TripolarGrid(arch;
                        size = (360, 180, 60),
                        z,
                        halo = (5, 5, 4))

    bottom_height = regrid_bathymetry(grid;
                                      minimum_depth = 10,
                                      major_basins = 2,
                                      interpolation_passes = 10)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height);
                                active_cells_map = true)

    free_surface = SplitExplicitFreeSurface(grid; substeps=70)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            free_surface,
                            closure,
                            kwargs...)
end
