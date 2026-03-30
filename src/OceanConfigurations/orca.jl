"""
    orca_ocean(arch = CPU(); zstar=false, kwargs...)

Construct an ocean `Simulation` on an ORCA grid (NEMO eORCA mesh) with
realistic bathymetry loaded from NEMO mesh-mask files.
"""
function orca_ocean(arch = CPU();
                    zstar = false,
                    momentum_advection = WENOVectorInvariant(order=5),
                    tracer_advection = WENO(order=5),
                    closure = default_one_degree_closure(),
                    kwargs...)

    z = vertical_coordinate(; zstar)

    grid = ORCAGrid(arch;
                    dataset = ORCA1(),
                    Nz = 60,
                    z,
                    halo = (4, 4, 4),
                    with_bathymetry = true,
                    active_cells_map = true)

    free_surface = SplitExplicitFreeSurface(grid; substeps=70)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            free_surface,
                            closure,
                            kwargs...)
end
