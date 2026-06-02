function default_sixth_degree_closure()
    catke = default_ocean_closure()
    vertical_diffusivity = VerticalScalarDiffusivity(ν=1e-5, κ=henyey_diffusivity)
    return (catke, vertical_diffusivity)
end

"""
    sixth_degree_tripolar_ocean(arch = CPU(); zstar=true, kwargs...)

Construct an ocean `Simulation` on a 1/6° (2160×1080) `TripolarGrid` with
realistic bathymetry. Designed for distributed multi-GPU runs
(e.g. `Distributed(GPU(), partition=Partition(2, 2))` on 4 GPUs).
"""
function sixth_degree_tripolar_ocean(arch = CPU();
                                     zstar = true,
                                     Nz = 60,
                                     depth = 6000,
                                     momentum_advection = WENOVectorInvariant(),
                                     tracer_advection = WENO(order=7),
                                     closure = default_sixth_degree_closure(),
                                     halo = (7, 7, 7),
                                     minimum_depth = 20,
                                     interpolation_passes = 40,
                                     substeps = 300,
                                     z = nothing,
                                     additional_surface_fluxes = nothing,
                                     kwargs...)

    if isnothing(z)
        z = vertical_coordinate(; Nz, depth, zstar)
    end

    grid = TripolarGrid(arch;
                        size = (2160, 1080, Nz),
                        z,
                        halo)

    bottom_height = regrid_bathymetry(grid;
                                      minimum_depth,
                                      major_basins = 1,
                                      interpolation_passes)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height);
                                active_cells_map = true)

    free_surface = SplitExplicitFreeSurface(grid; substeps)

    asf = resolve_surface_fluxes(additional_surface_fluxes, arch, grid)
    flux_kw = isnothing(asf) ? (;) : (; additional_surface_fluxes = asf)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            free_surface,
                            closure,
                            flux_kw...,
                            kwargs...)
end
