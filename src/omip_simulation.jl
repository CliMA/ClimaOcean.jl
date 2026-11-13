#####
##### The OMIP simulation
#####

@inline biharmonic_viscosity(i, j, k, grid, ℓx, ℓy, ℓz, clock, fields, λ) = Az(i, j, k, grid, ℓx, ℓy, ℓz)^2 / λ

horizontal_viscosity(::Nothing) = nothing
horizontal_viscosity(timescale) = HorizontalScalarBiharmonicDiffusivity(ν = biharmonic_viscosity, discrete_form = true, parameters = timescale)

function ocean_closure(ocean)
    background_mixing = VerticalScalarDiffusivity(κ = ocean.background_diffusivity, ν = ocean.background_viscosity)
    closures = (ocean.vertical_mixing, ocean.eddy_closure, horizontal_viscosity(ocean.biharmonic_timescale),
                background_mixing)
    return filter(!isnothing, closures)
end

# Search far enough for the river mouths to reach the ocean through about 3ᵒ of coastline at any resolution
function maximum_river_search_radius(grid)
    Nx, Ny, _ = size(grid)
    return max(5, ceil(Int, 3 / ((360 / Nx + 180 / Ny) / 2)))
end

"""
    omip_grid(configuration, arch)

The NEMO eORCA mesh of `configuration`, keeping only the largest connected ocean basin.
"""
omip_grid(configuration, arch) =
    ORCAGrid(arch; dataset = configuration.grid.dataset, Nz = configuration.grid.Nz,
             z = configuration.grid.vertical_coordinate, halo = (8, 8, 8), major_basins = 1,
             minimum_depth = configuration.grid.minimum_depth)

"""
    omip_ocean(grid, configuration; river_routing = nothing)

The ocean simulation of `configuration` on `grid`, starting from its initial condition converted to Conservative
Temperature and Absolute Salinity, and the conservative surface salinity restoring it is forced with.
"""
function omip_ocean(grid, configuration; river_routing = nothing)
    ocean = configuration.ocean

    chlorophyll = FieldTimeSeries(Metadata(:chlorophyll; dataset = ocean.chlorophyll_dataset,
                                           dates = (DateTime(2000, 1, 1), DateTime(2000, 12, 1))), grid;
                                  time_indices_in_memory = 12)

    top_cell_thickness = Δzᶜᶜᶜ(1, 1, size(grid, 3), host_grid(grid))

    salinity_restoring = DatasetRestoring(Metadata(:salinity; dataset = ocean.salinity_restoring_dataset),
                                          architecture(grid);
                                          rate = ocean.salinity_restoring_piston_velocity / top_cell_thickness,
                                          time_indices_in_memory = 12)

    convert_to_absolute_salinity!(salinity_restoring.field_time_series)
    salinity_flux = ConservativeSurfaceFluxRestoring(SurfaceFluxRestoring(salinity_restoring), grid)

    simulation = ocean_simulation(grid;
                                  momentum_advection = ocean.momentum_advection,
                                  tracer_advection = (T = ocean.tracer_advection, S = ocean.tracer_advection),
                                  free_surface = SplitExplicitFreeSurface(grid; substeps = ocean.barotropic_substeps),
                                  closure = ocean_closure(ocean),
                                  reference_density = ocean.reference_density,
                                  radiative_forcing = TwoColorRadiation(grid; chlorophyll),
                                  river_routing,
                                  additional_surface_fluxes = (; S = salinity_flux))

    T = CenterField(grid)
    S = CenterField(grid)
    set!(T, Metadatum(:temperature; dataset = ocean.initial_condition_dataset))
    set!(S, Metadatum(:salinity; dataset = ocean.initial_condition_dataset))
    convert_to_teos10!(T, S)
    flatten_polar_cap!(T, S, configuration.grid.polar_cap_latitude)
    set!(simulation.model; T, S)

    return simulation, salinity_flux
end

"""
    omip_sea_ice(grid, ocean, configuration)

The sea ice simulation of `configuration` on `grid`, coupled to `ocean` and starting from its initial condition.
"""
function omip_sea_ice(grid, ocean, configuration)
    sea_ice = configuration.sea_ice
    dynamics = sea_ice_dynamics(grid, ocean; sea_ice_ocean_drag_reference_depth = sea_ice.ocean_drag_reference_depth)
    itd_shape = UniformThicknessDistribution(categories = sea_ice.thickness_categories)

    simulation = sea_ice_simulation(grid, ocean; dynamics,
                                    advection = sea_ice.advection,
                                    internal_heat_flux = ConductiveFlux(; conductivity = sea_ice.conductivity, itd_shape),
                                    snow_thermodynamics = default_snow_thermodynamics(grid; itd_shape))

    set!(simulation.model, sea_ice.initial_condition)

    return simulation
end

"""
    omip_simulation(configuration = OneDegree();
                    arch = GPU(),
                    dates = (DateTime(1958, 1, 1), DateTime(2018, 1, 1)),
                    stop_time = Inf,
                    forcing_directory = NumericalEarth.DataWrangling.JRA55.download_JRA55_cache,
                    backend_size = 50,
                    diagnostics = true,
                    output_directory = ".",
                    kwargs...)

A coupled ocean--sea ice simulation in the configuration of the Ocean Model Intercomparison Project (OMIP-2)
protocol of [Griffies et al. (2016)](@cite griffies2016omip) and [Tsujino et al. (2020)](@cite tsujino2020evaluation),
forced by the JRA55-do reanalysis over `dates`. Every physical and numerical choice comes from the
[`OMIPConfiguration`](@ref).

Keyword arguments
=================

- `arch`: the architecture, possibly `Distributed`.
- `dates`: the `(start_date, end_date)` of the forcing.
- `stop_time`: the stop time of the simulation.
- `forcing_directory`: where the JRA55-do forcing is stored.
- `backend_size`: number of forcing records held in memory.
- `diagnostics`: whether to attach the output of [`add_omip_diagnostics!`](@ref).
- `output_directory`: where the output and checkpoints are written.
- remaining keyword arguments go to [`add_omip_diagnostics!`](@ref).
"""
function omip_simulation(configuration = OneDegree();
                         arch = GPU(),
                         dates = (DateTime(1958, 1, 1), DateTime(2018, 1, 1)),
                         stop_time = Inf,
                         forcing_directory = NumericalEarth.DataWrangling.JRA55.download_JRA55_cache,
                         backend_size = 50,
                         diagnostics = true,
                         output_directory = ".",
                         diagnostics_kw...)

    start_date, end_date = dates
    forcing = configuration.forcing
    interfaces = configuration.interfaces
    forcing_kw = (; dataset = forcing.dataset, start_date, end_date, dir = forcing_directory,
                    time_indices_in_memory = backend_size)

    grid = omip_grid(configuration, arch)

    land = JRA55PrescribedLand(grid; forcing_kw...,
                               maximum_search_radius = maximum_river_search_radius(grid),
                               spread_radius = forcing.river_spread_radius,
                               maximum_spread_cells = forcing.maximum_river_spread_cells)

    ocean, salinity_flux = omip_ocean(grid, configuration; river_routing = land.river_routing)
    sea_ice = omip_sea_ice(grid, ocean, configuration)

    atmosphere = JRA55PrescribedAtmosphere(arch; forcing_kw...)

    sea_ice_albedo = SeaIceAlbedo(sea_ice.model.ice_thickness,
                                  sea_ice.model.snow_thickness,
                                  sea_ice.model.snow_thermodynamics.top_surface_temperature)

    radiation = JRA55PrescribedRadiation(arch; forcing_kw...,
                                         ocean_surface = SurfaceRadiationProperties(interfaces.ocean_albedo, 1),
                                         sea_ice_surface = SurfaceRadiationProperties(sea_ice_albedo, 1))

    sea_ice_ocean_heat_flux = ThreeEquationHeatFlux(sea_ice;
                                                    heat_transfer_coefficient = interfaces.sea_ice_ocean_heat_transfer_coefficient,
                                                    friction_velocity = MomentumBasedFrictionVelocity())

    coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, land, radiation, sea_ice_ocean_heat_flux,
                                     interfaces.atmosphere_ocean_fluxes,
                                     interfaces.atmosphere_sea_ice_fluxes,
                                     interfaces.atmosphere_ocean_interface_specific_humidity,
                                     interfaces.atmosphere_sea_ice_interface_specific_humidity)

    time_step = configuration.time_step
    simulation = Simulation(coupled_model; Δt = time_step.initial_Δt, stop_time)

    add_callback!(simulation, ConservativeSurfaceFluxRestoringCallback(salinity_flux, ocean.model),
                  IterationInterval(1); name = :salinity_restoring)

    configuration.ocean.total_water_normalization &&
        add_callback!(simulation, TotalWaterNormalization(coupled_model), IterationInterval(1); name = :total_water)

    add_callback!(simulation, time_step, IterationInterval(1); name = :time_step)
    add_callback!(simulation, OMIPProgress(), IterationInterval(100); name = :progress)

    diagnostics && add_omip_diagnostics!(simulation; output_directory, prefix = configuration.name, diagnostics_kw...)

    return simulation
end

#####
##### Progress
#####

mutable struct OMIPProgress
    wall_time :: UInt64
    iteration :: Int
end

OMIPProgress() = OMIPProgress(time_ns(), 0)

Oceananigans.prognostic_state(::OMIPProgress) = nothing

function (progress::OMIPProgress)(simulation)
    ocean = simulation.model.ocean.model
    sea_ice = simulation.model.sea_ice.model
    u, v, w = ocean.velocities
    T, S = ocean.tracers.T, ocean.tracers.S

    seconds_per_step = 1e-9 * (time_ns() - progress.wall_time) / max(1, iteration(simulation) - progress.iteration)

    @info @sprintf("%s, iteration %d, Δt %s, max|u| (%.2f, %.2f, %.2e) m s⁻¹, T ∈ [%.1f, %.1f] ᵒC, S ∈ [%.1f, %.1f] g kg⁻¹, max hⁱ %.1f m, %s per step",
                   prettytime(simulation), iteration(simulation), prettytime(simulation.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   minimum(T), maximum(T), minimum(S), maximum(S),
                   maximum(sea_ice.ice_thickness), prettytime(seconds_per_step))

    progress.wall_time = time_ns()
    progress.iteration = iteration(simulation)

    return nothing
end
