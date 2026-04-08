using Printf

"""
    omip_simulation(config::Symbol = :half_degree; kwargs...)

Create a fully coupled ocean--sea-ice--atmosphere OMIP simulation.

The single `config` argument selects the grid configuration:

- `:half_degree` -- 720x360 tripolar grid with EN4Monthly restoring
- `:orca`        -- NEMO eORCA mesh with WOAMonthly restoring

Returns a `Simulation` wrapping an `OceanSeaIceModel`.

# Keyword arguments

- `arch`: architecture (`CPU()` or `GPU()`). Default: `CPU()`.
- `Nz::Int`: number of vertical levels. Default: `100`.
- `depth`: maximum ocean depth in metres. Default: `5500`.
- `κ_skew`: skew diffusivity for the GM closure. Default: `500`.
- `κ_symmetric`: symmetric diffusivity for the GM closure. Default: `100`.
- `forcing_dir`: directory for JRA55 forcing data. Default: `"forcing_data"`.
- `restoring_dir`: directory for restoring / initial-condition climatology. Default: `"climatology"`.
- `restoring_rate`: surface salinity restoring piston velocity in m/day. Default: `1/6`.
- `start_date`: simulation start date. Default: `DateTime(1958, 1, 1)`.
- `end_date`: end date for forcing / restoring metadata. Default: `DateTime(1958, 12, 30)`.
- `restart`: path to a JLD2 checkpoint file, or `nothing`. Default: `nothing`.
- `Δt`: simulation time step. Default: `30minutes`.
- `stop_time`: simulation stop time. Default: `Inf`.
- `diagnostics::Bool`: whether to attach OMIP diagnostics. Default: `true`.
- `surface_averaging_interval`: averaging window for surface diagnostics. Default: `15days`.
- `field_averaging_interval`: averaging window for 3-D field diagnostics. Default: `15days`.
- `checkpoint_interval`: interval between checkpoint writes. Default: `90days`.
- `output_dir`: directory for output files. Default: `"."`.
- `filename_prefix`: prefix for output file names. Default: `string(config)`.
- `file_splitting_interval`: time interval for splitting output files. Default: `365days` (1 year).
- `transports`: whether to include transport diagnostics. Default: `false`.
"""
function omip_simulation(config::Symbol = :half_degree;
                         arch = CPU(),
                         Nz = 100,
                         depth = 5500,
                         κ_skew = 500,
                         κ_symmetric = 100,
                         forcing_dir = "forcing_data",
                         restoring_dir = "climatology",
                         restoring_rate = 1/6,   # m/day
                         start_date = DateTime(1958, 1, 1),
                         end_date = DateTime(1958, 12, 30),
                         restart = nothing,
                         Δt = 30minutes,
                         stop_time = Inf,
                         diagnostics = true,
                         surface_averaging_interval = 365days / 24,
                         field_averaging_interval = 365days / 24,
                         checkpoint_interval = 90days,
                         output_dir = ".",
                         filename_prefix = string(config),
                         file_splitting_interval = 365days,
                         transports = false)

    cfg = Val(config)

    # 1. Build the ocean simulation (grid + physics + salinity restoring + ICs)
    ocean = _build_ocean(cfg, arch;
                         Nz, depth, κ_skew, κ_symmetric,
                         restoring_dir, restoring_rate,
                         start_date, end_date)

    grid = ocean.model.grid

    # 2. Build the sea-ice simulation and set initial conditions
    sea_ice = _build_sea_ice(cfg, grid, ocean; restoring_dir)

    # 3. Build the prescribed atmosphere
    atmosphere, radiation = omip_atmosphere(arch;
                                           forcing_dir,
                                           start_date,
                                           end_date)

    # 4. Couple into an OceanSeaIceModel
    coupled = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

    # 5. Optionally load a restart checkpoint
    if !isnothing(restart)
        _load_restart!(coupled, atmosphere, ocean, sea_ice, restart)
    end

    # 6. Wrap in a Simulation
    simulation = Simulation(coupled; Δt, stop_time)

    # 7. Attach progress callback
    wall_time = Ref(time_ns())
    add_callback!(simulation, _omip_progress_callback(wall_time), IterationInterval(10))

    # 8. Optionally attach diagnostics
    if diagnostics
        add_omip_diagnostics!(simulation;
                              surface_averaging_interval,
                              field_averaging_interval,
                              checkpoint_interval,
                              output_dir,
                              filename_prefix,
                              file_splitting_interval)

        if transports
            add_omip_transport_diagnostics!(simulation;
                                            output_dir,
                                            filename_prefix,
                                            field_averaging_interval)
        end
    end

    return simulation
end

#####
##### Ocean builders
#####

function _build_ocean(::Val{:half_degree}, arch;
                      Nz, depth, κ_skew, κ_symmetric,
                      restoring_dir, restoring_rate,
                      start_date, end_date)

    tracer_advection = WENO(order=7, minimum_buffer_upwind_order=3)

    # Build salinity restoring forcing before ocean construction,
    # so it can be passed directly to the ocean simulation constructor.
    FS = _salinity_restoring_forcing(arch, Nz, depth, EN4Monthly();
                                     restoring_dir, restoring_rate,
                                     start_date, end_date)

    ocean = half_degree_tripolar_ocean(arch;
                                       Nz, depth,
                                       κ_skew, κ_symmetric,
                                       tracer_advection,
                                       forcing = (; S = FS))

    # Initial conditions from EN4Monthly
    set!(ocean.model,
         T = Metadatum(:temperature; dir=restoring_dir, dataset=EN4Monthly(), date=start_date),
         S = Metadatum(:salinity;    dir=restoring_dir, dataset=EN4Monthly(), date=start_date))

    return ocean
end

function _build_ocean(::Val{:orca}, arch;
                      Nz, depth, κ_skew, κ_symmetric,
                      restoring_dir, restoring_rate,
                      start_date, end_date)

    FS = _salinity_restoring_forcing(arch, Nz, depth, WOAMonthly();
                                     restoring_dir, restoring_rate,
                                     start_date, end_date)

    ocean = orca_ocean(arch;
                       Nz, depth,
                       κ_skew, κ_symmetric,
                       forcing = (; S = FS))

    # Initial conditions from WOAMonthly
    set!(ocean.model,
         T = Metadatum(:temperature; dir=restoring_dir, dataset=WOAMonthly(), date=start_date),
         S = Metadatum(:salinity;    dir=restoring_dir, dataset=WOAMonthly(), date=start_date))

    return ocean
end

#####
##### Salinity restoring (shared by both configurations)
#####

"""
    _salinity_restoring_forcing(arch, Nz, depth, dataset; kw...)

Build a `DatasetRestoring` forcing term for surface salinity restoring.

The restoring is applied only to the top grid cell. The `restoring_rate` (m/day)
is converted to a volumetric rate by dividing by the top-cell thickness.
"""
function _salinity_restoring_forcing(arch, Nz, depth, dataset;
                                     restoring_dir,
                                     restoring_rate,
                                     start_date,
                                     end_date)

    # Compute the vertical coordinate to determine surface-cell geometry.
    # Uses the same ExponentialDiscretization as OceanConfigurations.vertical_coordinate.
    z = ExponentialDiscretization(Nz, -depth, 0)

    # z(Nz+1) is the top face (should be 0), z(Nz) is the face below it.
    z_surface  = z(Nz) - 1  # match reference: znode(Nz, grid, Face()) - 1
    Δz_surface = z(Nz + 1) - z(Nz)

    rate = restoring_rate / (Δz_surface * days)

    surface_mask(x, y, z, t) = z >= z_surface

    Smetadata = Metadata(:salinity;
                         dir = restoring_dir,
                         dataset,
                         start_date,
                         end_date)

    FS = DatasetRestoring(Smetadata, arch;
                          rate,
                          mask = surface_mask,
                          time_indices_in_memory = 12)

    return FS
end

#####
##### Sea-ice builders
#####

function _build_sea_ice(::Val{:half_degree}, grid, ocean; restoring_dir)
    sea_ice = half_degree_tripolar_sea_ice(ocean)

    set!(sea_ice.model,
         h = Metadatum(:sea_ice_thickness;     dir=restoring_dir, dataset=ECCO4Monthly()),
         ℵ = Metadatum(:sea_ice_concentration; dir=restoring_dir, dataset=ECCO4Monthly()))

    return sea_ice
end

function _build_sea_ice(::Val{:orca}, grid, ocean; restoring_dir)
    sea_ice = orca_sea_ice(ocean)

    set!(sea_ice.model,
         h = Metadatum(:sea_ice_thickness;     dir=restoring_dir, dataset=ECCO4Monthly()),
         ℵ = Metadatum(:sea_ice_concentration; dir=restoring_dir, dataset=ECCO4Monthly()))

    return sea_ice
end

#####
##### Restart loading
#####

"""
    _load_restart!(coupled, atmosphere, ocean, sea_ice, checkpoint_path)

Restore model state from a JLD2 checkpoint file. Synchronizes clocks across
all components (ocean, sea ice, atmosphere, coupled model).
"""
function _load_restart!(coupled, atmosphere, ocean, sea_ice, checkpoint_path)
    ocean_model  = ocean.model
    seaice_model = sea_ice.model

    file = jldopen(checkpoint_path)

    try
        uo = file["uo"]
        vo = file["vo"]
        wo = file["wo"]
        To = file["To"]
        So = file["So"]
        ηo = file["ηo"]
        ui = file["ui"]
        vi = file["vi"]
        hi = file["hi"]
        ℵi = file["ℵi"]

        clock = file["clock"]

        set!(ocean_model, u=uo, v=vo, w=wo, T=To, S=So, η=ηo)

        # Restore TKE if present
        eo = try
            file["eo"]
        catch
            nothing
        end

        if !isnothing(eo)
            try
                set!(ocean_model, e=eo)
            catch
                @warn "Could not restore TKE (e) from checkpoint"
            end
        end

        set!(seaice_model, h=hi, ℵ=ℵi)
        set!(seaice_model.velocities.u, ui)
        set!(seaice_model.velocities.v, vi)

        # Synchronize all clocks to the checkpoint clock
        _synch_clock!(ocean_model.clock, clock)
        _synch_clock!(seaice_model.clock, clock)

        Oceananigans.initialize!(ocean_model)

        # Synchronize the atmosphere and coupled-model clocks
        _synch_clock!(atmosphere, ocean_model)
        _synch_clock!(coupled, ocean_model)

        time_step!(atmosphere, 0)
        atmosphere.clock.iteration -= 1

        @info "Restarted from checkpoint: $(checkpoint_path)"
        @info "  Ocean clock:       $(ocean_model.clock)"
        @info "  Sea-ice clock:     $(seaice_model.clock)"
        @info "  Atmosphere clock:  $(atmosphere.clock)"
        @info "  Coupled clock:     $(coupled.clock)"
    finally
        close(file)
    end

    return nothing
end

# Clock synchronization utilities

function _synch_clock!(clock1::Oceananigans.TimeSteppers.Clock, clock2)
    clock1.time      = clock2.time
    clock1.iteration = clock2.iteration
    clock1.last_Δt   = clock2.last_Δt
    return nothing
end

_synch_clock!(model1, model2) = _synch_clock!(model1.clock, model2.clock)

#####
##### Progress callback
#####

function _omip_progress_callback(wall_time)
    function progress(sim)
        sea_ice = sim.model.sea_ice
        ocean   = sim.model.ocean

        hmax = maximum(sea_ice.model.ice_thickness)
        ℵmax = maximum(sea_ice.model.ice_concentration)
        Tmax = maximum(ocean.model.tracers.T)
        Tmin = minimum(ocean.model.tracers.T)
        Smax = maximum(ocean.model.tracers.S)
        Smin = minimum(ocean.model.tracers.S)
        umax = maximum(ocean.model.velocities.u)
        vmax = maximum(ocean.model.velocities.v)
        wmax = maximum(ocean.model.velocities.w)

        step_time = 1e-9 * (time_ns() - wall_time[])

        msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ",
                        prettytime(sim), iteration(sim), prettytime(sim.Δt))
        msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
        msg3 = @sprintf("extrema(T, S): (%.2f, %.2f) C, (%.2f, %.2f) psu ",
                        Tmin, Tmax, Smin, Smax)
        msg4 = @sprintf("maximum(u): (%.2e, %.2e, %.2e) m/s, ", umax, vmax, wmax)
        msg5 = @sprintf("wall time: %s", prettytime(step_time))

        @info msg1 * msg2 * msg3 * msg4 * msg5

        wall_time[] = time_ns()

        return nothing
    end

    return progress
end
