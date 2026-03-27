function default_sixth_degree_closure()
    catke = default_ocean_closure()
    horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true, parameters=60days)
    vertical_diffusivity = VerticalScalarDiffusivity(ν=1e-5, κ=henyey_diffusivity)
    return (catke, horizontal_viscosity, vertical_diffusivity)
end

"""
    sixth_degree_tripolar_ocean(arch = CPU(); zstar=false, kwargs...)

Construct an ocean `Simulation` on a 1/6° (2160×1080) `TripolarGrid` with
realistic bathymetry. Designed for distributed multi-GPU runs
(e.g. `Distributed(GPU(), partition=Partition(2, 2))` on 4 GPUs).
"""
function sixth_degree_tripolar_ocean(arch = CPU();
                                     zstar = false,
                                     momentum_advection = WENOVectorInvariant(order=5),
                                     tracer_advection = WENO(order=7),
                                     closure = default_sixth_degree_closure(),
                                     kwargs...)

    z = vertical_coordinate(; zstar)

    grid = TripolarGrid(arch;
                        size = (2160, 1080, 60),
                        z,
                        halo = (7, 7, 7))

    bottom_height = regrid_bathymetry(grid;
                                      minimum_depth = 20,
                                      major_basins = 1,
                                      interpolation_passes = 40)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height);
                                active_cells_map = true)

    free_surface = SplitExplicitFreeSurface(grid; substeps=300)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            free_surface,
                            closure,
                            kwargs...)
end
