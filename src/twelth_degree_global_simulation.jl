using Oceananigans.Utils 
using Oceananigans.MultiRegion: reconstruct_global_field, multi_region_object_from_array
using Oceananigans.TurbulenceClosures: HorizontalDivergenceFormulation
using Oceananigans.Advection: VelocityStencil

"""
    twelth_degree_near_global_simulation(architecture = GPU(); kwargs...)

Return an Oceananigans.Simulation of Earth's ocean at 1/12 degree resolution.
The data for initial conditions, boundary conditions and bathymetry _MUST_ be locally available (momentarily too large for online storage)
"""
function twelth_degree_near_global_simulation(architecture = GPU();
        size                                         = (4320, 1800, 48),
        number_of_partions                           = 4,
        boundary_layer_turbulence_closure            = RiBasedVerticalDiffusivity(),
        background_vertical_diffusivity              = 1e-5,
        background_vertical_viscosity                = 1e-4,
        horizontal_viscosity                         = leith_viscosity(HorizontalDivergenceFormulation(), C_vort = 2.0, C_div = 3.0),
        surface_temperature_relaxation_time_scale    = 7days,
        surface_salinity_relaxation_time_scale       = 7days,
        bottom_drag_coefficient                      = 3e-3,
        reference_density                            = 1029.0,  
        time_step                                    = 2minutes,
        stop_iteration                               = Inf,
        start_time                                   = 345days,
        stop_time                                    = Inf,
        equation_of_state                            = TEOS10EquationOfState(; reference_density),
        tracers                                      = [:T, :S],
        initial_conditions                           = "../data/evolved-initial-conditions-165days.jld2",
        bathymetry_path                              = "../data/bathymetry-ad-hoc.jld2",
        boundary_conditions_path                     = "../data/boundary_conditions_twelth_degree.jld2"
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
    τˣ =  - jldopen(boundary_conditions_path)["τˣ"] ./ reference_density
    τʸ =  - jldopen(boundary_conditions_path)["τʸ"] ./ reference_density
    T★ =    jldopen(boundary_conditions_path)["Tₛ"] 
    S★ =    jldopen(boundary_conditions_path)["Sₛ"] 
    F★ =    zeros(Base.size(S★)...)
    Q★ =    zeros(Base.size(T★)...)
    close(boundary_conditions_path)
    @info "... read boundary conditions (" * prettytime(1e-9 * (time_ns() - start)) * ")"

    # Stretched faces from ECCO Version 4 (49 levels in the vertical)
    z_faces = VerticalGrids.z_49_levels_10_to_400_meter_spacing

    # A spherical domain
    underlying_grid = LatitudeLongitudeGrid(architecture; size,
                                            longitude = (-180, 180),
                                            latitude = (-75, 75),
                                            halo = (5, 5, 5),
                                            z = z_faces)

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))

    underlying_mrg = MultiRegionGrid(underlying_grid, partition = XPartition(number_of_partions), devices = number_of_partitions);
    mrg            = MultiRegionGrid(grid,            partition = XPartition(number_of_partions), devices = number_of_partitions);
    
    @info "Created $grid"

    # Convert boundary conditions arrays to multi GPU
    τˣ = multi_region_object_from_array(τˣ, mrg);
    τʸ = multi_region_object_from_array(τʸ, mrg);
    target_sea_surface_temperature = T★ = multi_region_object_from_array(T★, mrg);
    target_sea_surface_salinity    = S★ = multi_region_object_from_array(S★, mrg);
    surface_temperature_flux       = Q★ = multi_region_object_from_array(Q★, mrg);
    surface_salt_flux              = F★ = multi_region_object_from_array(F★, mrg);

    #####
    ##### Physics and model setup
    #####

    vitd = VerticallyImplicitTimeDiscretization()

    vertical_viscosity   = VerticalScalarDiffusivity(vitd, ν=background_vertical_viscosity, κ=background_vertical_diffusivity)

    closures = Any[horizontal_viscosity, boundary_layer_turbulence_closure, vertical_viscosity]

    boundary_layer_turbulence_closure isa CATKEVerticalDiffusivity &&
        push!(tracers, :e)

    # TODO: do this internally in model constructor
    closures = tuple(closures...)

    #####
    ##### Boundary conditions / time-dependent fluxes 
    #####

    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form=true, parameters = bottom_drag_coefficient)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form=true, parameters = bottom_drag_coefficient)

    no_slip_bc = ValueBoundaryCondition(0)

    u_immersed_bc = ImmersedBoundaryCondition(bottom = drag_u,
                                              west = no_slip_bc,
                                              east = no_slip_bc,
                                              south = no_slip_bc,
                                              north = no_slip_bc)

    v_immersed_bc = ImmersedBoundaryCondition(bottom = drag_v,
                                              west = no_slip_bc,
                                              east = no_slip_bc,
                                              south = no_slip_bc,
                                              north = no_slip_bc)

    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)

    Nmonths = 12 # number of months in the forcing file
    u_wind_stress_parameters = (; τ=τˣ, Nmonths)
    v_wind_stress_parameters = (; τ=τʸ, Nmonths)
    u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form=true, parameters=u_wind_stress_parameters)
    v_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form=true, parameters=v_wind_stress_parameters)

    Δz_top = @allowscalar Δzᵃᵃᶜ(1, 1, grid.Nz, grid.underlying_grid)

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

    buoyancy = SeawaterBuoyancy(; equation_of_state)
    coriolis = HydrostaticSphericalCoriolis(scheme = WetCellEnstrophyConservingScheme())
    free_surface = ImplicitFreeSurface()

    model = HydrostaticFreeSurfaceModel(; grid = mrg, free_surface, coriolis, buoyancy, tracers,
                                          momentum_advection = WENO(vector_invariant = VelocityStencil()),
                                          closure = closures,
                                          boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
                                          tracer_advection = WENO(underlying_mrg))

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
    
        u = reconstruct_global_field(sim.model.velocities.u)
        w = reconstruct_global_field(sim.model.velocities.w)
    
        intw  = Array(interior(w))
        max_w = findmax(intw)
    
        mw = max_w[1]
        iw = max_w[2]
    
        @info @sprintf("Time: % 12s, iteration: %d, max(|u|): %.2e ms⁻¹, wmax: %.2e , loc: (%d, %d, %d), wall time: %s", 
                        prettytime(sim.model.clock.time),
                        sim.model.clock.iteration, maximum(abs, u), mw, iw[1], iw[2], iw[3], 
                        prettytime(wall_time))
    
        start_time[1] = time_ns()
    
        return nothing
    end
    
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    return simulation
end
