"""
    orca_ocean(arch = CPU(); zstar=true, Nz=60, depth=6000, kwargs...)

Construct an ocean `Simulation` on an ORCA grid (NEMO eORCA mesh) with
realistic bathymetry loaded from NEMO mesh-mask files.
"""
function orca_ocean(arch = CPU();
                    zstar = true,
                    Nz = 60,
                    depth = 6000,
                    momentum_advection = WENOVectorInvariant(order=5),
                    tracer_advection = WENO(order=5),
                    κ_skew = 500,
                    κ_symmetric = 200,
                    biharmonic_timescale = 15days,
                    background_κ = henyey_diffusivity,
                    background_ν = 1e-5,
                    halo = (4, 4, 4),
                    substeps = 70,
                    kwargs...)

    z = vertical_coordinate(; Nz, depth, zstar)

    closure = default_one_degree_closure(; κ_skew, κ_symmetric, biharmonic_timescale, background_κ, background_ν)

    grid = ORCAGrid(arch;
                    dataset = ORCA1(),
                    Nz,
                    z,
                    halo,
                    with_bathymetry = true,
                    active_cells_map = true)

    free_surface = SplitExplicitFreeSurface(grid; substeps)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            free_surface,
                            closure,
                            kwargs...)
end
