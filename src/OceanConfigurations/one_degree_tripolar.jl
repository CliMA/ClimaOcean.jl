function default_one_degree_closure(; κ_skew=500,
                                     κ_symmetric=200,
                                     biharmonic_timescale=15days,
                                     background_κ=henyey_diffusivity,
                                     background_ν=1e-5)
    catke = default_ocean_closure()
    eddy  = IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric)
    horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true, parameters=biharmonic_timescale)
    vertical_diffusivity = VerticalScalarDiffusivity(ν=background_ν, κ=background_κ)
    return (catke, eddy, horizontal_viscosity, vertical_diffusivity)
end

"""
    one_degree_tripolar_ocean(arch = CPU(); zstar=false, Nz=60, depth=6000, kwargs...)

Construct an ocean `Simulation` on a one-degree (360×180) `TripolarGrid` with
realistic bathymetry and closures suitable for coarser resolution
(CATKE + Gent-McWilliams + biharmonic viscosity).
"""
function one_degree_tripolar_ocean(arch = CPU();
                                   zstar = false,
                                   Nz = 60,
                                   depth = 6000,
                                   momentum_advection = WENOVectorInvariant(order=5),
                                   tracer_advection = WENO(order=5),
                                   κ_skew = 500,
                                   κ_symmetric = 200,
                                   biharmonic_timescale = 15days,
                                   background_κ = henyey_diffusivity,
                                   background_ν = 1e-5,
                                   halo = (5, 5, 4),
                                   minimum_depth = 10,
                                   interpolation_passes = 10,
                                   substeps = 70,
                                   kwargs...)

    z = vertical_coordinate(; Nz, depth, zstar)

    closure = default_one_degree_closure(; κ_skew, κ_symmetric, biharmonic_timescale, background_κ, background_ν)

    grid = TripolarGrid(arch;
                        size = (360, 180, Nz),
                        z,
                        halo)

    bottom_height = regrid_bathymetry(grid;
                                      minimum_depth,
                                      major_basins = 2,
                                      interpolation_passes)

    grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height);
                                active_cells_map = true)

    free_surface = SplitExplicitFreeSurface(grid; substeps)

    return ocean_simulation(grid;
                            momentum_advection,
                            tracer_advection,
                            free_surface,
                            closure,
                            kwargs...)
end
