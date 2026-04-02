using Oceananigans.Units
using JLD2

"""
    add_omip_diagnostics!(simulation; kwargs...)

Attach OMIP-protocol output writers to a coupled ocean--sea-ice simulation.

Creates three output writers:
1. **Surface diagnostics** (`{prefix}_surface.jld2`): 2-D fields averaged over
   `surface_averaging_interval` — SST, SSS, SSH, surface velocity, squared
   fields for variance, mixed-layer depth, wind stress, heat/freshwater fluxes,
   and sea-ice state variables.
2. **3-D field diagnostics** (`{prefix}_fields.jld2`): Full 3-D temperature,
   salinity, velocity, and (if available) TKE, averaged over
   `field_averaging_interval`.
3. **Checkpointer** (`{prefix}_checkpoint`): JLD2 checkpoint at
   `checkpoint_interval`.

Optionally, transport diagnostics (MOC, heat transports) can be enabled with
`transports = true`.

# Keyword arguments

- `surface_averaging_interval`: averaging window for surface output. Default: `365days / 24`.
- `field_averaging_interval`: averaging window for 3-D output. Default: `365days / 24`.
- `checkpoint_interval`: interval between checkpoints. Default: `90days`.
- `output_dir`: directory for all output files. Default: `"."`.
- `filename_prefix`: prefix for output filenames. Default: `"omip"`.
- `file_splitting_interval`: time interval for splitting output files. Default: `365days` (1 year).
"""
function add_omip_diagnostics!(simulation;
                               surface_averaging_interval = 365days / 24,
                               field_averaging_interval = 365days / 24,
                               checkpoint_interval = 90days,
                               output_dir = ".",
                               filename_prefix = "omip",
                               file_splitting_interval = 365days)

    file_splitting = TimeInterval(file_splitting_interval)

    model    = simulation.model          # OceanSeaIceModel
    ocean    = model.ocean               # ocean Simulation
    sea_ice  = model.sea_ice             # sea-ice Simulation
    grid     = ocean.model.grid
    Nz       = size(grid, 3)

    # ─── Ocean fields ───────────────────────────────────────────────────
    T, S = ocean.model.tracers.T, ocean.model.tracers.S
    u, v, w = ocean.model.velocities
    η = ocean.model.free_surface.displacement

    # ─── Flux fields from coupled-model interfaces ──────────────────────
    τx = model.interfaces.net_fluxes.ocean.u       # zonal wind stress
    τy = model.interfaces.net_fluxes.ocean.v       # meridional wind stress
    JT = model.interfaces.net_fluxes.ocean.T       # net heat flux
    Js = model.interfaces.net_fluxes.ocean.S       # net freshwater flux
    Qc = model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
    Qv = model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat

    # ─── Sea-ice fields ─────────────────────────────────────────────────
    hi = sea_ice.model.ice_thickness
    ℵi = sea_ice.model.ice_concentration
    ui, vi = sea_ice.model.velocities

    # Sea-ice top surface temperature (may not be available for all models)
    sitemptop = try
        sea_ice.model.ice_thermodynamics.top_surface_temperature
    catch
        nothing
    end

    # ─── Mixed-layer depth ──────────────────────────────────────────────
    mld = MixedLayerDepthField(ocean.model.buoyancy, grid, ocean.model.tracers)

    # ===================================================================
    # 1. Surface diagnostics
    # ===================================================================
    surface_indices = (:, :, Nz)

    tos = Field(T; indices = surface_indices)       # SST
    sos = Field(S; indices = surface_indices)       # SSS

    uo_surface = Field(u; indices = surface_indices)
    vo_surface = Field(v; indices = surface_indices)

    # Squared surface fields for computing variance
    tossq = Field(T * T; indices = surface_indices)
    sossq = Field(S * S; indices = surface_indices)
    zossq = η * η

    surface_outputs = Dict{Symbol, Any}(
        :tos        => tos,
        :sos        => sos,
        :zos        => η,
        :uo_surface => uo_surface,
        :vo_surface => vo_surface,
        :tossq      => tossq,
        :sossq      => sossq,
        :zossq      => zossq,
        :mlotst     => mld,
        :tauuo      => τx,
        :tauvo      => τy,
        :hfds       => JT,
        :wfo        => Js,
        :hfss       => Qc,
        :hfls       => Qv,
        :siconc     => ℵi,
        :sithick    => hi,
        :siu        => ui,
        :siv        => vi,
    )

    if !isnothing(sitemptop)
        surface_outputs[:sitemptop] = sitemptop
    end

    simulation.output_writers[:surface] = JLD2Writer(ocean.model, surface_outputs;
        schedule = AveragedTimeInterval(surface_averaging_interval),
        filename = joinpath(output_dir, filename_prefix * "_surface"),
        file_splitting,
        overwrite_existing = true)

    # ===================================================================
    # 2. 3-D field diagnostics
    # ===================================================================
    field_outputs = Dict{Symbol, Any}(
        :thetao => T,
        :so     => S,
        :uo     => u,
        :vo     => v,
        :wo     => w,
    )

    # Include TKE if the turbulence closure provides it
    if haskey(ocean.model.tracers, :e)
        field_outputs[:tke] = ocean.model.tracers.e
    end

    simulation.output_writers[:fields] = JLD2Writer(ocean.model, field_outputs;
        schedule = AveragedTimeInterval(field_averaging_interval),
        filename = joinpath(output_dir, filename_prefix * "_fields"),
        file_splitting,
        overwrite_existing = true)

    # ===================================================================
    # 3. Checkpointer
    # ===================================================================
    simulation.output_writers[:checkpointer] = Checkpointer(simulation.model;
        schedule = TimeInterval(checkpoint_interval),
        prefix   = joinpath(output_dir, filename_prefix * "_checkpoint"),
        cleanup  = false,
        verbose  = true)

    # ===================================================================
    # 4. Scalar / zonal-mean diagnostics (callback-based)
    # ===================================================================
    scalar_writes_per_file = max(1, floor(Int, file_splitting_interval / field_averaging_interval))
    scalar_cb = OMIPScalarCallback(grid, T, S,
                                   joinpath(output_dir, filename_prefix * "_scalars.jld2");
                                   max_writes_per_file = scalar_writes_per_file)

    add_callback!(simulation, scalar_cb, TimeInterval(field_averaging_interval))

    @info "OMIP diagnostics attached:" *
          " surface ($(length(surface_outputs)) fields, every $(prettytime(surface_averaging_interval)))," *
          " 3-D ($(length(field_outputs)) fields, every $(prettytime(field_averaging_interval)))," *
          " scalars/zonal means (every $(prettytime(field_averaging_interval)))," *
          " checkpointer (every $(prettytime(checkpoint_interval)))"

    return nothing
end

#####
##### Scalar / zonal-mean callback
#####

"""
    OMIPScalarCallback

Callback that computes zonal means of T and S (latitude x depth) and
global scalar diagnostics, appending results to a JLD2 file at each invocation.

Saved fields:
- `T_zonal_mean`, `S_zonal_mean`: zonal-mean profiles (latitude x depth)
- `tosga`, `sosga`: global-mean SST and SSS
- `thetaoga`, `soga`: global-mean T and S (full volume)
- `latitude`, `depth`: coordinate arrays (saved once on first call)
"""
mutable struct OMIPScalarCallback{G, FT, FS}
    grid :: G
    T_field :: FT
    S_field :: FS
    filepath :: String
    max_writes_per_file :: Int
    writes_in_current_file :: Int
    part :: Int
    initialized :: Bool
end

OMIPScalarCallback(grid, T, S, filepath; max_writes_per_file=24) =
    OMIPScalarCallback(grid, T, S, filepath, max_writes_per_file, 0, 1, false)

function (cb::OMIPScalarCallback)(simulation)
    grid = cb.grid
    T = cb.T_field
    S = cb.S_field
    Nz = size(grid, 3)
    t = Oceananigans.time(simulation)
    iter = iteration(simulation)

    # Zonal means (reuses infrastructure from Diagnostics module)
    T̄, S̄, φ, z = compute_zonal_averages(grid, T, S)

    # Global-mean surface fields
    SST = Array(interior(T, :, :, Nz))
    SSS = Array(interior(S, :, :, Nz))
    tosga = mean(filter(isfinite, SST))
    sosga = mean(filter(isfinite, SSS))

    # Global-mean full-volume fields
    T_all = Array(interior(T))
    S_all = Array(interior(S))
    thetaoga = mean(filter(isfinite, T_all))
    soga     = mean(filter(isfinite, S_all))

    # Rotate file when max writes reached
    if cb.writes_in_current_file >= cb.max_writes_per_file
        cb.part += 1
        cb.writes_in_current_file = 0
        cb.initialized = false  # new file needs coordinates
    end

    # Write to JLD2
    key = "timeseries/iteration$(iter)"

    jldopen(_scalar_filepath(cb), "a+") do file
        if !cb.initialized
            file["latitude"] = Float64.(φ)
            file["depth"]    = Float64.(z)
            cb.initialized = true
        end

        file["$(key)/time"]          = t
        file["$(key)/T_zonal_mean"]  = Float32.(T̄)
        file["$(key)/S_zonal_mean"]  = Float32.(S̄)
        file["$(key)/tosga"]         = Float32(tosga)
        file["$(key)/sosga"]         = Float32(sosga)
        file["$(key)/thetaoga"]      = Float32(thetaoga)
        file["$(key)/soga"]          = Float32(soga)
    end

    cb.writes_in_current_file += 1

    return nothing
end

function _scalar_filepath(cb::OMIPScalarCallback)
    if cb.part == 1
        return cb.filepath
    else
        base, ext = splitext(cb.filepath)
        return "$(base)_part$(cb.part)$(ext)"
    end
end
