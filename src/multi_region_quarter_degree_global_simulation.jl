using Oceananigans.TurbulenceClosures: HorizontalDivergenceFormulation
using Oceananigans.Advection: VelocityStencil
using Oceananigans.MultiRegion: multi_region_object_from_array

"""
    quarter_degree_near_global_simulation(architecture = GPU(); kwargs...)

Return an Oceananigans.Simulation of Earth's ocean at 1/4 degree resolution.
"""
function quarter_degree_near_global_simulation(architecture = GPU();
        size                                         = (1440, 600, 48),
        boundary_layer_turbulence_closure            = RiBasedVerticalDiffusivity(),
        background_vertical_diffusivity              = 1e-5,
        background_vertical_viscosity                = 1e-4,
        horizontal_viscosity                         = geometric_viscosity(HorizontalDivergenceFormulation(), 5days),
        surface_temperature_relaxation_time_scale    = 7days,
        surface_salinity_relaxation_time_scale       = 7days,
        bottom_drag_coefficient                      = 3e-3,
        reference_density                            = 1029.0,
        reference_heat_capacity                      = 3991.0,
        reference_salinity                           = 34.0,
        time_step                                    = 6minutes,
        stop_iteration                               = Inf,
        start_time                                   = 345days,
        stop_time                                    = Inf,
        equation_of_state                            = TEOS10EquationOfState(; reference_density),
        tracers                                      = [:T, :S],
        initial_conditions                           = datadep"near_global_quarter_degree/initial_conditions.jld2",
        bathymetry_path                              = datadep"near_global_quarter_degree/bathymetry-1440x600.jld2",
        temp_surface_boundary_conditions_path        = datadep"near_global_quarter_degree/temp-1440x600-latitude-75.jld2",
        salt_surface_boundary_conditions_path        = datadep"near_global_quarter_degree/salt-1440x600-latitude-75.jld2",
        u_stress_surface_boundary_conditions_path    = datadep"near_global_quarter_degree/tau_x-1440x600-latitude-75.jld2",
        v_stress_surface_boundary_conditions_path    = datadep"near_global_quarter_degree/tau_y-1440x600-latitude-75.jld2",
)

    bathymetry_file = jldopen(bathymetry_path)
    bathymetry = bathymetry_file["bathymetry"]
    close(bathymetry_file)

    @info "Reading initial conditions..."; start=time_ns()
    initial_conditions_file = jldopen(initial_conditions)
    T_init = initial_conditions_file["T"]
    S_init = initial_conditions_file["S"]
    close(initial_conditions_file)
    @info "... read initial conditions (" * prettytime(1e-9 * (time_ns() - start)) * ")"

    # Files contain 12 arrays of monthly-averaged data from 1992
    @info "Reading boundary conditions..."; start=time_ns()
    # Files contain 1 year (1992) of 12 monthly averages
    τˣ =  - jldopen(u_stress_surface_boundary_conditions_path)["field"] ./ reference_density
    τʸ =  - jldopen(v_stress_surface_boundary_conditions_path)["field"] ./ reference_density
    T★ =    jldopen(temp_surface_boundary_conditions_path)["field"] 
    S★ =    jldopen(salt_surface_boundary_conditions_path)["field"] 
    F★ =    zeros(Base.size(S★)...)
    Q★ =    zeros(Base.size(T★)...)
    
    @info "... read boundary conditions (" * prettytime(1e-9 * (time_ns() - start)) * ")"

    # Stretched faces from ECCO Version 4 (49 levels in the vertical)
    z_faces = VerticalGrids.z_49_levels_10_to_400_meter_spacing

    # A spherical domain
    underlying_grid = LatitudeLongitudeGrid(architecture; size,
                                            longitude = (-180, 180),
                                            latitude = (-75, 75),
                                            halo = (5, 5, 5),
                                            z = z_faces)
    
    Δz_top = @allowscalar Δzᵃᵃᶜ(1, 1, underlying_grid.Nz, underlying_grid)

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))
    underlying_grid = MultiRegionGrid(underlying_grid, partition = XPartition(4), devices = 4)
    grid            = MultiRegionGrid(grid,            partition = XPartition(4), devices = 4)

    # Convert boundary conditions arrays to GPU
    τˣ = multi_region_object_from_array(τˣ, grid)
    τʸ = multi_region_object_from_array(τʸ, grid)
    target_sea_surface_temperature = T★ = multi_region_object_from_array(architecture, T★)
    target_sea_surface_salinity    = S★ = multi_region_object_from_array(architecture, S★)
    surface_temperature_flux       = Q★ = multi_region_object_from_array(architecture, Q★)
    surface_salt_flux              = F★ = multi_region_object_from_array(architecture, F★)

    T_init = multi_region_object_from_array(T_init, grid)
    S_init = multi_region_object_from_array(S_init, grid)

    @info "Created $grid"

    #####
    ##### Physics and model setup
    #####

    vitd = VerticallyImplicitTimeDiscretization()

    vertical_viscosity   = VerticalScalarDiffusivity(vitd, ν=background_vertical_viscosity, κ=background_vertical_diffusivity)

    closures = Any[boundary_layer_turbulence_closure, vertical_viscosity]

    boundary_layer_turbulence_closure isa CATKEVerticalDiffusivity &&
        push!(tracers, :e)

    # TODO: do this internally in model constructor
    closures = tuple(closures...)

    #####
    ##### Boundary conditions / time-dependent fluxes 
    #####

    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = bottom_drag_coefficient)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = bottom_drag_coefficient)

    u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u)
    v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v)

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

    Nmonths = 12 # number of months in the forcing file
    u_wind_stress_parameters = (; τ=τˣ, Nmonths)
    v_wind_stress_parameters = (; τ=τʸ, Nmonths)
    u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form=true, parameters=u_wind_stress_parameters)
    v_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form=true, parameters=v_wind_stress_parameters)

    T_relaxation_parameters = (; λ = Δz_top / surface_temperature_relaxation_time_scale,
                                 Nmonths,
                                 T★ = target_sea_surface_temperature,
                                 Q★ = surface_temperature_flux)

    S_relaxation_parameters = (; λ = Δz_top / surface_salinity_relaxation_time_scale,
                                 Nmonths,
                                 S★ = target_sea_surface_salinity,
                                 F★ = surface_salt_flux)

    T_surface_relaxation_bc = FluxBoundaryCondition(surface_temperature_relaxation,
                                                    discrete_form = true,
                                                    parameters = T_relaxation_parameters)

    S_surface_relaxation_bc = FluxBoundaryCondition(surface_salinity_relaxation,
                                                    discrete_form = true,
                                                    parameters = S_relaxation_parameters)

    u_bcs = FieldBoundaryConditions(top = u_wind_stress_bc,
                                    bottom = u_bottom_drag_bc,
                                    immersed = u_immersed_bc)

    v_bcs = FieldBoundaryConditions(top = v_wind_stress_bc,
                                    bottom = v_bottom_drag_bc,
                                    immersed = v_immersed_bc)

    T_bcs = FieldBoundaryConditions(top = T_surface_relaxation_bc)
    S_bcs = FieldBoundaryConditions(top = S_surface_relaxation_bc)

    buoyancy     = SeawaterBuoyancy(; equation_of_state)
    coriolis     = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme())
    free_surface = ImplicitFreeSurface()

    model = HydrostaticFreeSurfaceModel(; grid, free_surface, coriolis, buoyancy, tracers,
                                          momentum_advection = VectorInvariant(vorticity_scheme  = WENO(),
                                                                               divergence_schem  = WENO(),
                                                                               vertical_scheme   = WENO(underlying_grid)),
                                          closure = closures,
                                          boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
                                          tracer_advection = WENO(underlying_grid))

    @info "... built $model."
    @info "Model building time: " * prettytime(1e-9 * (time_ns() - start))

    #####
    ##### Initial condition:
    #####

    set!(model, T=T_init, S=S_init)

    # Because MITgcm forcing starts at Jan 15 (?)
    model.clock.time = start_time

    simulation = Simulation(model; Δt=time_step, stop_iteration, stop_time)

    start_time = [time_ns()]

    function progress(sim)
        wall_time = (time_ns() - start_time[1]) * 1e-9

        u = sim.model.velocities.u
        w = sim.model.velocities.w

        intw  = Array(interior(w))
        max_w = findmax(intw)

        mw = max_w[1]
        iw = max_w[2]

        msg1 = @sprintf("Time: % 12s, iteration: %d, ", prettytime(sim), iteration(sim))
        msg2 = @sprintf("max(|u|): %.2e ms⁻¹, wmax: %.2e, loc: (%d, %d, %d), ",
                        maximum(abs, u), mw, iw[1], iw[2], iw[3])
        msg3 = @sprintf("wall time: %s", prettytime(wall_time))

        @info msg1 * msg2 * msg3

        start_time[1] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    return simulation
end
