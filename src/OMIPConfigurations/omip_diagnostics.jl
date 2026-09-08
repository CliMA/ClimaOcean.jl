
using JLD2
using Oceananigans.AbstractOperations: KernelFunctionOperation, Integral
using Oceananigans.Operators: Axᶠᶜᶜ, Ayᶜᶠᶜ, ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, δzᵃᵃᶜ, Δz⁻¹ᶜᶠᶜ
using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity, getclosure, κ_ϵSyᶜᶠᶠ

# Volumetric face fluxes [m³/s] for offline transport diagnostics
# (e.g. AMOC streamfunction). Saved as 3-D fields on the velocity faces
# so a downstream zonal sum + cumulative-z gives the streamfunction
# directly without re-deriving Δx/Δz on the visualizer side.
@inline zonal_volume_flux(i, j, k, grid, u)      = @inbounds u[i, j, k] * Axᶠᶜᶜ(i, j, k, grid)
@inline meridional_volume_flux(i, j, k, grid, v) = @inbounds v[i, j, k] * Ayᶜᶠᶜ(i, j, k, grid)

# GM eddy-induced meridional velocity v★ = -δz(κ ϵ Sʸ) / Δz, for the `DiffusiveFormulation`, where GM
# enters as a skew tracer flux and v★ is never materialized. Built from the same helpers Oceananigans
# uses for the `AdvectiveFormulation`, so the two formulations give bitwise-identical diagnostics.
@inline function skew_meridional_velocity(i, j, k, grid, clock, closure, buoyancy, model_fields)
    closure = getclosure(i, j, closure)
    κ = closure.κ_skew
    limiter = closure.slope_limiter
    return - δzᵃᵃᶜ(i, j, k, grid, κ_ϵSyᶜᶠᶠ, clock, limiter, κ, buoyancy, model_fields) * Δz⁻¹ᶜᶠᶜ(i, j, k, grid)
end

"""
    bolus_meridional_velocity(model)

Return the GM eddy-induced meridional velocity — the field the `AdvectiveFormulation` already advects
with, or an operation rebuilding it under the `DiffusiveFormulation`. `nothing` without a skew
diffusivity, so the eddy-resolving configurations write no identically-zero field to disk.
"""
function bolus_meridional_velocity(model)
    closures = model.closure isa Tuple ? model.closure : (model.closure,)
    K        = model.closure isa Tuple ? model.closure_fields : (model.closure_fields,)

    n = findfirst(c -> c isa IsopycnalSkewSymmetricDiffusivity && !isnothing(c.κ_skew), closures)
    isnothing(n) && return nothing

    hasproperty(K[n], :v) && return K[n].v

    return KernelFunctionOperation{Center, Face, Center}(skew_meridional_velocity, model.grid,
                                                         model.clock, closures[n], model.buoyancy,
                                                         fields(model))
end

# Square-then-interpolate helpers for the centred velocity-squared
# diagnostics (`uosq`, `vosq`, `kega`). Same pattern as `ϕ²` /
# `τᶜᶜᶜ` in `friction_velocity.jl` and the KE loader in
# `scripts/visualize/cache.jl` — written as `KernelFunctionOperation`s
# so the surrounding `Average(...)` reduction (which `compute!`s the
# operand at every averaging timestep) runs a single fused kernel
# rather than the multi-kernel AbstractOperation tree the previous
# `@at((Center, Center, Center), u*u)` form built.
@inline ψ²(i, j, k, grid, ψ) = @inbounds ψ[i, j, k]^2

@inline uu_at_ccc(i, j, k, grid, u) = ℑxᶜᵃᵃ(i, j, k, grid, ψ², u)
@inline vv_at_ccc(i, j, k, grid, v) = ℑyᵃᶜᵃ(i, j, k, grid, ψ², v)

@inline ke_at_ccc(i, j, k, grid, u, v) =
    (ℑxᶜᵃᵃ(i, j, k, grid, ψ², u) + ℑyᵃᶜᵃ(i, j, k, grid, ψ², v)) / 2

# Unit operand so `Integral` returns ∫dV = ocean volume (z-star-live, immersed cells
# carry zero metric). Closes the freshwater budget: volume drift = net water added.
@inline one_at_ccc(i, j, k, grid) = one(grid)

"""
    add_omip_diagnostics!(simulation; kwargs...)

Attach OMIP-protocol output writers to a coupled ocean--sea-ice
simulation built by [`omip_simulation`](@ref).

Creates four output writers:

1. **Surface diagnostics** (`<prefix>_surface.nc`): 2-D fields averaged
   over `surface_averaging_interval` -- SST, SSS, SSH, surface velocities,
   squared fields for variance, mixed-layer depth, wind stress,
   heat/freshwater fluxes, and sea-ice state variables.
2. **3-D field diagnostics** (`<prefix>_fields.nc`): full 3-D temperature,
   salinity, velocity, buoyancy, and (when present) TKE, averaged over
   `field_averaging_interval`.
3. **Averages** (`<prefix>_averages.nc`): global means of T, S, buoyancy
   and horizontal-mean (dims=(1,2)) depth profiles of the same, on the
   same `field_averaging_interval` schedule.
4. **Checkpointer** (`<prefix>_checkpoint`): JLD2 checkpoint of the
   coupled model at `checkpoint_interval`. Use `run!(sim; pickup=true)`
   to restart from the latest checkpoint.

# Keyword arguments

- `surface_averaging_interval`: averaging window for surface output. Default: `5days`.
- `field_averaging_interval`: averaging window for 3-D / averages output. Default: `15days`.
- `checkpoint_interval`: interval between checkpoints. Default: `90days`.
- `output_dir`: directory for all output files. Default: `"."`.
- `filename_prefix`: prefix for output filenames. Default: `"omip"`.
- `file_splitting_interval`: time interval for splitting output files. Default: `360days`.
"""
function add_omip_diagnostics!(simulation;
                               field_mean_interval = 5days,
                               surface_averaging_interval = 5days,
                               field_averaging_interval = 15days,
                               checkpoint_interval = 720days,
                               output_dir = ".",
                               filename_prefix = "omip",
                               file_splitting_interval = 360days)

    model    = simulation.model
    ocean    = model.ocean
    sea_ice  = model.sea_ice
    grid     = ocean.model.grid
    Nz       = size(grid, 3)

    T, S = ocean.model.tracers.T, ocean.model.tracers.S
    u, v, w = ocean.model.velocities
    η = ocean.model.free_surface.displacement

    τx = model.interfaces.net_fluxes.ocean.u
    τy = model.interfaces.net_fluxes.ocean.v
    JT = model.interfaces.net_fluxes.ocean.T
    Js = model.interfaces.net_fluxes.ocean.S
    Qc = model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
    Qv = model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat

    JTf  = NumericalEarth.Diagnostics.frazil_temperature_flux(model)
    JTn  = NumericalEarth.Diagnostics.net_ocean_temperature_flux(model)
    JTio = NumericalEarth.Diagnostics.sea_ice_ocean_temperature_flux(model)
    JTao = NumericalEarth.Diagnostics.atmosphere_ocean_temperature_flux(model)
    JSn  = NumericalEarth.Diagnostics.net_ocean_salinity_flux(model)
    JSio = NumericalEarth.Diagnostics.sea_ice_ocean_salinity_flux(model)

    hi = sea_ice.model.ice_thickness
    ℵi = sea_ice.model.ice_concentration
    ui, vi = sea_ice.model.velocities

    sitemptop = try
        sea_ice.model.ice_thermodynamics.top_surface_temperature
    catch
        nothing
    end

    mld = MixedLayerDepthField(ocean.model.buoyancy, grid, ocean.model.tracers)

    # Surface diagnostics
    surface_indices = (:, :, Nz)

    tos = view(T, :, :, Nz)
    sos = view(S, :, :, Nz)
    uo_surface = view(u, :, :, Nz)
    vo_surface = view(v, :, :, Nz)

    tossq = tos * tos
    sossq = sos * sos
    zossq = Field(η * η)

    surface_outputs = Dict{Symbol, Any}(
        :tos      => tos,
        :sos      => sos,
        :zos      => η,
        :uos      => uo_surface,
        :vos      => vo_surface,
        :tossq    => tossq,
        :sossq    => sossq,
        :zossq    => zossq,
        :mlotst   => mld,
        :tauuo    => τx,
        :tauvo    => τy,
        :hfds     => JT,
        :wfo      => Js,
        :hfss     => Qc,
        :hfls     => Qv,
        :siconc   => ℵi,
        :sithick  => hi,
        :siu      => ui,
        :siv      => vi,
        :JTf      => JTf,
        :JTn      => JTn,
        :JTio     => JTio,
        :JTao     => JTao,
        :JSn      => JSn,
        :JSio     => JSio
    )

    if !isnothing(sitemptop)
        surface_outputs[:sitemptop] = sitemptop
    end

    hs = sea_ice.model.snow_thickness
    if !isnothing(hs)
        surface_outputs[:sisnthick] = hs
    end

    simulation.output_writers[:surface] = JLD2Writer(ocean.model, surface_outputs;
                                                     schedule = AveragedTimeInterval(surface_averaging_interval),
                                                     dir = output_dir,
                                                     filename = filename_prefix * "_surface",
                                                     file_splitting = TimeInterval(file_splitting_interval),
                                                     overwrite_existing = true,
                                                     jld2_kw = Dict(:compress => ZstdFilter()))

    # 3-D fields (including buoyancy)
    bop = Oceananigans.Models.buoyancy_operation(ocean.model)

    uvol = KernelFunctionOperation{Face,   Center, Center}(zonal_volume_flux,      grid, u)
    vvol = KernelFunctionOperation{Center, Face,   Center}(meridional_volume_flux, grid, v)

    uosq = KernelFunctionOperation{Center, Center, Center}(uu_at_ccc, grid, u)
    vosq = KernelFunctionOperation{Center, Center, Center}(vv_at_ccc, grid, v)

    field_outputs = Dict{Symbol, Any}(
        :to   => T,
        :so   => S,
        :uo   => u,
        :vo   => v,
        :wo   => w,
        :bo   => bop,
        :uosq => uosq,
        :vosq => vosq,
        :uvol => uvol,
        :vvol => vvol,
    )

    if haskey(ocean.model.tracers, :e)
        field_outputs[:tke] = ocean.model.tracers.e
    end

    # Same ×Aʸ wrapper as `vvol`, so the two add directly offline into the residual streamfunction.
    # Column-integrating this flux telescopes to zero, which is why the bolus moves the overturning
    # but leaves full-depth strait transports untouched.
    vgm = bolus_meridional_velocity(ocean.model)
    if !isnothing(vgm)
        field_outputs[:vvolgm] = KernelFunctionOperation{Center, Face, Center}(meridional_volume_flux, grid, vgm)
    end

    # Field-valued (flow-aware) eddy coefficients: write what the closure actually used, so a run
    # with KSKEW=nemo/cesm can be audited against the offline climatological estimates.
    closures = ocean.model.closure isa Tuple ? ocean.model.closure : (ocean.model.closure,)
    n = findfirst(c -> c isa IsopycnalSkewSymmetricDiffusivity, closures)
    if !isnothing(n)
        closures[n].κ_skew isa Field      && (field_outputs[:aei] = closures[n].κ_skew)
        closures[n].κ_symmetric isa Field && (field_outputs[:aht] = closures[n].κ_symmetric)
    end

    simulation.output_writers[:fields] = JLD2Writer(ocean.model, field_outputs;
                                                    schedule = AveragedTimeInterval(field_averaging_interval),
                                                    dir = output_dir,
                                                    filename = filename_prefix * "_fields",
                                                    file_splitting = TimeInterval(file_splitting_interval),
                                                    overwrite_existing = true,
                                                    jld2_kw = Dict(:compress => ZstdFilter()))

    # `:zosga` (area-mean free-surface displacement) is a Boussinesq mass-conservation check;
    # the horizontal metric is fixed, so the 2-D `Average` is exact. Depth profiles are per-level
    # area-means (also exact). Global volume-means, in contrast, are stored as exact `Integral`s and
    # divided by `:voco` in the visualizer, because `Average` freezes its volume denominator at
    # construction — wrong on a z-star grid. `:soco`/`:hoco` double as salt/heat content (×ρ/1000
    # for salt kg, ×ρ·cp for heat J); `:sivol` closes the ocean+ice salt budget via the ice reservoir.
    # `:sivol`/`:snvol` also close the *water* budget: freshwater normalization corrects only the
    # atmospheric part of the surface flux, so `:voco` alone is expected to move with the sea-ice
    # exchange. The conserved quantity is ocean + ice + snow, weighted by their densities.
    vol_op = KernelFunctionOperation{Center, Center, Center}(one_at_ccc, grid)
    average_outputs = Dict{Symbol, Any}(
        :zosga => Average(η),
        :to_h  => Average(T,   dims=(1, 2)),
        :so_h  => Average(S,   dims=(1, 2)),
        :bo_h  => Average(bop, dims=(1, 2)),
        :voco  => Integral(vol_op),
        :soco  => Integral(S),
        :hoco  => Integral(T),
        :sivol => Integral(hi * ℵi),
    )

    if !isnothing(hs)
        average_outputs[:snvol] = Integral(hs * ℵi)
    end

    simulation.output_writers[:averages] = JLD2Writer(ocean.model, average_outputs;
                                                      schedule = AveragedTimeInterval(field_mean_interval),
                                                      dir = output_dir,
                                                      filename = filename_prefix * "_averages",
                                                      file_splitting = TimeInterval(file_splitting_interval),
                                                      overwrite_existing = true)

    # Checkpointer (drives `run!(sim; pickup=true)`)
    simulation.output_writers[:checkpointer] = Checkpointer(simulation.model;
                                                            schedule = TimeInterval(checkpoint_interval),
                                                            prefix   = joinpath(output_dir, filename_prefix * "_checkpoint"),
                                                            cleanup  = false,
                                                            verbose  = true)

    @info "OMIP diagnostics attached:" *
          " surface ($(length(surface_outputs)) fields, every $(prettytime(surface_averaging_interval)))," *
          " 3-D ($(length(field_outputs)) fields, every $(prettytime(field_averaging_interval)))," *
          " averages ($(length(average_outputs)) fields, every $(prettytime(field_averaging_interval)))," *
          " checkpointer (every $(prettytime(checkpoint_interval)))"

    return nothing
end
