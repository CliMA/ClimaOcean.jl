# OMIP Configurations Pipeline — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a `ClimaOcean.OMIPConfigurations` module that provides a turnkey `omip_simulation(:half_degree)` / `omip_simulation(:orca)` API with standardized OMIP diagnostics, producing NetCDF output and JLD2 checkpoints.

**Architecture:** The module composes on top of the existing `OceanConfigurations` and `SeaIceConfigurations` factories, passing OMIP-specific physical parameters. Diagnostics are attached via a standalone `add_omip_diagnostics!` function that can also be used independently on any simulation. Output is split into surface (2D), fields (3D), scalars (zonal/global means), and checkpoints (JLD2).

**Tech Stack:** Oceananigans.jl, NumericalEarth.jl, ClimaSeaIce.jl, NCDatasets.jl (for NetCDF output)

---

## Background & Reference

### OMIP Protocol (Griffies et al. 2016)

The OMIP diagnostic protocol prescribes Priority 1 outputs (all monthly by default, we use 15-day default):

**Surface 2D fields:** SST (tos), SSS (sos), SSH (zos), bottom T/S (tob/sob), MLD (mlotst), squared fields (tossq, sossq, zossq), surface velocity, barotropic streamfunction

**3D fields:** temperature (thetao), salinity (so), velocity (uo, vo, wo), buoyancy frequency (obvfsq), cell thickness

**Boundary fluxes:** net heat flux (hfds) and components (rsntds, rlntds, hfls, hfss), net water flux (wfo) and components (pr, evs, friver), wind stress (tauuo, tauvo), sea ice salt flux

**Scalars:** global mean T, S, SST, SSS, ocean mass, ocean volume

**Transport diagnostics (opt-in):** MOC streamfunction, meridional heat/salt transport, strait mass transports

### Current codebase state

- `src/OceanConfigurations/half_degree_tripolar.jl:16-46` — factory function, accepts `closure` kwarg as full closure tuple
- `src/OceanConfigurations/orca.jl:7-32` — same pattern, uses `default_one_degree_closure()`
- `src/OceanConfigurations/OceanConfigurations.jl:21-27` — shared `νhb` and `henyey_diffusivity` functions
- `src/Diagnostics/` — has `MixedLayerDepthField`, `compute_report_fields`, zonal average infrastructure
- `NumericalEarth.jl/experiments/omips/omip_defaults.jl` — the monolithic OMIP factory to be replaced (255 lines, commented-out diagnostics at lines 104-187 show the flux access patterns)

### Key flux access patterns (from omip_defaults.jl:173-181)
```julia
τx = model.interfaces.net_fluxes.ocean.u          # wind stress x
τy = model.interfaces.net_fluxes.ocean.v          # wind stress y
JT = model.interfaces.net_fluxes.ocean.T          # net heat flux
Js = model.interfaces.net_fluxes.ocean.S          # net salt flux
Qc = model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
Qv = model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat
Qi = model.interfaces.sea_ice_ocean_interface.fluxes.interface_heat
Ji = model.interfaces.sea_ice_ocean_interface.fluxes.salt
```

---

## Task 1: Refactor `OceanConfigurations` to accept physical parameters

**Files:**
- Modify: `src/OceanConfigurations/OceanConfigurations.jl`
- Modify: `src/OceanConfigurations/half_degree_tripolar.jl`
- Modify: `src/OceanConfigurations/orca.jl`
- Modify: `src/OceanConfigurations/one_degree_tripolar.jl`

**Goal:** Replace the `closure` kwarg with physical parameter kwargs. The factories build closures internally from parameters. This lets OMIP (and any user) customize physics without constructing closure objects.

### Step 1: Refactor `half_degree_tripolar.jl`

Replace the current content with:

```julia
"""
    default_half_degree_closure(; κ_skew=500, κ_symmetric=200,
                                  biharmonic_timescale=40days,
                                  background_κ=henyey_diffusivity,
                                  background_ν=1e-5)

Build the standard closure tuple for half-degree ocean simulations.
"""
function default_half_degree_closure(; κ_skew = 500,
                                       κ_symmetric = 200,
                                       biharmonic_timescale = 40days,
                                       background_κ = henyey_diffusivity,
                                       background_ν = 1e-5)

    catke = default_ocean_closure()
    eddy  = IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric)
    horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true,
                                                                  parameters=biharmonic_timescale)
    vertical_diffusivity = VerticalScalarDiffusivity(ν=background_ν, κ=background_κ)
    return (catke, eddy, horizontal_viscosity, vertical_diffusivity)
end

"""
    half_degree_tripolar_ocean(arch = CPU(); Nz=60, depth=6000,
                               κ_skew=500, κ_symmetric=200,
                               biharmonic_timescale=40days,
                               zstar=false, kwargs...)

Construct an ocean `Simulation` on a half-degree (720×360) `TripolarGrid`.
Physics are configured via physical parameters, not closure objects.
"""
function half_degree_tripolar_ocean(arch = CPU();
                                    Nz = 60,
                                    depth = 6000,
                                    zstar = false,
                                    κ_skew = 500,
                                    κ_symmetric = 200,
                                    biharmonic_timescale = 40days,
                                    background_κ = henyey_diffusivity,
                                    background_ν = 1e-5,
                                    momentum_advection = WENOVectorInvariant(order=5),
                                    tracer_advection = WENO(order=7),
                                    halo = (7, 7, 7),
                                    minimum_depth = 20,
                                    interpolation_passes = 25,
                                    substeps = 150,
                                    kwargs...)

    closure = default_half_degree_closure(; κ_skew, κ_symmetric,
                                            biharmonic_timescale,
                                            background_κ, background_ν)

    z = vertical_coordinate(; Nz, depth, zstar)

    grid = TripolarGrid(arch; size=(720, 360, Nz), z, halo)

    bottom_height = regrid_bathymetry(grid;
                                       minimum_depth,
                                       major_basins = 1,
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
```

**Key changes from current code:**
- `closure` kwarg removed, replaced by `κ_skew`, `κ_symmetric`, `biharmonic_timescale`, `background_κ`, `background_ν`
- New kwargs: `Nz` (was hardcoded 60), `depth` (was hardcoded 6000), `halo`, `minimum_depth`, `interpolation_passes`, `substeps`
- `default_half_degree_closure` now accepts kwargs
- `vertical_coordinate` now accepts `Nz` and `depth` (see Step 3)

### Step 2: Refactor `orca.jl`

Same pattern — expose physical parameters, add `Nz`, `depth`:

```julia
function orca_ocean(arch = CPU();
                    Nz = 60,
                    depth = 6000,
                    zstar = false,
                    κ_skew = 500,
                    κ_symmetric = 200,
                    biharmonic_timescale = 15days,
                    background_κ = henyey_diffusivity,
                    background_ν = 1e-5,
                    momentum_advection = WENOVectorInvariant(order=5),
                    tracer_advection = WENO(order=5),
                    substeps = 70,
                    kwargs...)

    closure = default_one_degree_closure(; κ_skew, κ_symmetric,
                                           biharmonic_timescale,
                                           background_κ, background_ν)

    z = vertical_coordinate(; Nz, depth, zstar)

    grid = ORCAGrid(arch;
                    dataset = ORCA1(),
                    Nz, z,
                    halo = (4, 4, 4),
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
```

### Step 3: Refactor `one_degree_tripolar.jl`

Apply the same parameter-based pattern to `default_one_degree_closure`:

```julia
function default_one_degree_closure(; κ_skew = 500,
                                      κ_symmetric = 200,
                                      biharmonic_timescale = 15days,
                                      background_κ = henyey_diffusivity,
                                      background_ν = 1e-5)
    catke = default_ocean_closure()
    eddy  = IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric)
    horizontal_viscosity = HorizontalScalarBiharmonicDiffusivity(ν=νhb, discrete_form=true,
                                                                  parameters=biharmonic_timescale)
    vertical_diffusivity = VerticalScalarDiffusivity(ν=background_ν, κ=background_κ)
    return (catke, eddy, horizontal_viscosity, vertical_diffusivity)
end
```

And update `one_degree_tripolar_ocean` to expose `Nz`, `depth`, and physics params (same pattern as half-degree).

### Step 4: Update `vertical_coordinate` in `OceanConfigurations.jl`

```julia
function vertical_coordinate(; Nz=60, depth=6000, zstar=false)
    return ExponentialDiscretization(Nz, -depth, 0; mutable=zstar)
end
```

### Step 5: Verify existing examples still work

Run: `julia --project -e 'using ClimaOcean; half_degree_tripolar_ocean()'` (should work with all defaults unchanged)

### Step 6: Commit

```
feat(OceanConfigurations): expose physical parameters as kwargs

Replace closure kwarg with individual physical parameters (κ_skew,
κ_symmetric, biharmonic_timescale, etc.) so that OMIP and other
configurations can customize physics without building closure objects.
Also expose Nz and depth as kwargs (previously hardcoded).
```

---

## Task 2: Create `OMIPConfigurations` module skeleton

**Files:**
- Create: `src/OMIPConfigurations/OMIPConfigurations.jl`
- Create: `src/OMIPConfigurations/omip_simulation.jl`
- Create: `src/OMIPConfigurations/atmosphere.jl`
- Create: `src/OMIPConfigurations/diagnostics.jl`
- Create: `src/OMIPConfigurations/transport_diagnostics.jl`
- Modify: `src/ClimaOcean.jl`

### Step 1: Create `OMIPConfigurations.jl` module file

```julia
module OMIPConfigurations

using Oceananigans
using Oceananigans.Units
using Dates

using NumericalEarth.Oceans: ocean_simulation
using NumericalEarth.SeaIces: sea_ice_simulation
using NumericalEarth.Atmospheres: JRA55PrescribedAtmosphere
using NumericalEarth.EarthSystemModels: OceanSeaIceModel, Radiation
using NumericalEarth.DataWrangling: Metadatum, Metadata, DatasetRestoring,
                                    EN4Monthly, ECCO4Monthly, WOAMonthly,
                                    MultiYearJRA55, JRA55NetCDFBackend

using ..OceanConfigurations: half_degree_tripolar_ocean, orca_ocean
using ..SeaIceConfigurations: half_degree_tripolar_sea_ice, orca_sea_ice
using ..Diagnostics: MixedLayerDepthField

export omip_simulation, add_omip_diagnostics!

include("atmosphere.jl")
include("diagnostics.jl")
include("transport_diagnostics.jl")
include("omip_simulation.jl")

end # module
```

**Note:** The exact import paths from NumericalEarth may need adjustment — verify at implementation time by checking what `NumericalEarth` actually exports.

### Step 2: Wire into `src/ClimaOcean.jl`

Add after line 81 (after `include("SeaIceConfigurations/SeaIceConfigurations.jl")`):

```julia
include("OMIPConfigurations/OMIPConfigurations.jl")
```

Add to `using` block after line 86:

```julia
using .OMIPConfigurations
```

Add to exports (after line 14):

```julia
       omip_simulation,
       add_omip_diagnostics!,
```

### Step 3: Commit

```
feat: add OMIPConfigurations module skeleton
```

---

## Task 3: Implement `omip_simulation`

**Files:**
- Write: `src/OMIPConfigurations/omip_simulation.jl`
- Write: `src/OMIPConfigurations/atmosphere.jl`

### Step 1: Implement `atmosphere.jl`

Shared atmosphere setup used by all OMIP configs:

```julia
"""
    omip_atmosphere(arch; forcing_dir, start_date, end_date, backend_size=30)

Set up JRA55 prescribed atmosphere with river and iceberg forcing.
"""
function omip_atmosphere(arch;
                         forcing_dir,
                         start_date,
                         end_date,
                         backend_size = 30)

    dataset = MultiYearJRA55()
    backend = JRA55NetCDFBackend(backend_size)

    atmosphere = JRA55PrescribedAtmosphere(arch;
                                           dir = forcing_dir,
                                           dataset,
                                           backend,
                                           include_rivers_and_icebergs = true,
                                           start_date,
                                           end_date)

    radiation = Radiation()

    return atmosphere, radiation
end
```

### Step 2: Implement `omip_simulation.jl`

```julia
"""
    omip_simulation(config::Symbol; kwargs...)

Build a complete OMIP coupled ocean-sea-ice simulation.

`config` can be `:half_degree` or `:orca`.

Returns a `Simulation` wrapping an `OceanSeaIceModel`, with diagnostics
attached by default.

# Keyword Arguments

## Grid and physics
- `arch = CPU()`: architecture (CPU or GPU)
- `Nz = 100`: number of vertical levels
- `depth = 5500`: ocean depth in meters
- `κ_skew = 500`: skew diffusivity for GM parameterization
- `κ_symmetric = 100`: symmetric (Redi) diffusivity

## Forcing and restoring
- `forcing_dir = "forcing_data"`: path to JRA55 atmospheric forcing files
- `restoring_dir = "climatology"`: path to restoring climatology data
- `restoring_rate = 1/6`: piston velocity for surface salinity restoring (m/day)
- `start_date = DateTime(1958, 1, 1)`: simulation start date
- `end_date = DateTime(1958, 12, 30)`: end date for forcing/restoring data

## Restart
- `restart = nothing`: path to checkpoint JLD2 file for restart

## Time stepping
- `Δt = 30minutes`: simulation time step
- `stop_time = Inf`: simulation stop time

## Diagnostics
- `diagnostics = true`: attach OMIP diagnostics (pass `false` or `nothing` to skip)
- `surface_averaging_interval = 15days`: averaging window for 2D surface output
- `field_averaging_interval = 15days`: averaging window for 3D field output
- `checkpoint_interval = 90days`: interval for JLD2 checkpoints
- `output_dir = "."`: directory for output files
- `filename_prefix = string(config)`: prefix for output filenames
- `transports = false`: include transport diagnostics (MOC, heat transport, straits)
"""
function omip_simulation(config::Symbol;
                         # Architecture
                         arch = CPU(),
                         # Grid
                         Nz = 100,
                         depth = 5500,
                         # Physics (passed through to OceanConfigurations)
                         κ_skew = 500,
                         κ_symmetric = 100,
                         # Forcing & restoring
                         forcing_dir = "forcing_data",
                         restoring_dir = "climatology",
                         restoring_rate = 1/6, # m/day
                         start_date = DateTime(1958, 1, 1),
                         end_date = DateTime(1958, 12, 30),
                         # Restart
                         restart = nothing,
                         # Time stepping
                         Δt = 30minutes,
                         stop_time = Inf,
                         # Diagnostics
                         diagnostics = true,
                         surface_averaging_interval = 15days,
                         field_averaging_interval = 15days,
                         checkpoint_interval = 90days,
                         output_dir = ".",
                         filename_prefix = string(config),
                         transports = false)

    # --- Build ocean ---
    ocean = _build_ocean(Val(config), arch; Nz, depth, κ_skew, κ_symmetric,
                         restoring_dir, restoring_rate, start_date, end_date)

    # --- Build sea ice ---
    grid = ocean.model.grid
    sea_ice = _build_sea_ice(Val(config), grid, ocean; restoring_dir)

    # --- Build atmosphere ---
    atmosphere, radiation = omip_atmosphere(arch; forcing_dir, start_date, end_date)

    # --- Couple ---
    coupled = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

    # --- Restart ---
    if !isnothing(restart)
        _load_restart!(coupled, atmosphere, ocean, sea_ice, restart)
    end

    # --- Simulation ---
    simulation = Simulation(coupled; Δt, stop_time)

    # --- Diagnostics ---
    if diagnostics === true || diagnostics === :default
        add_omip_diagnostics!(simulation;
                              surface_averaging_interval,
                              field_averaging_interval,
                              checkpoint_interval,
                              output_dir,
                              filename_prefix,
                              transports)
    end

    return simulation
end

# ----- Half-degree dispatch -----

function _build_ocean(::Val{:half_degree}, arch;
                      Nz, depth, κ_skew, κ_symmetric,
                      restoring_dir, restoring_rate, start_date, end_date)

    ocean = half_degree_tripolar_ocean(arch;
                                       Nz, depth,
                                       κ_skew, κ_symmetric,
                                       tracer_advection = WENO(order=7, minimum_buffer_upwind_order=3))

    grid = ocean.model.grid

    # Surface salinity restoring to EN4
    Smetadata = Metadata(:salinity; dir=restoring_dir, dataset=EN4Monthly(),
                         start_date, end_date)

    Nz_grid = size(grid, 3)
    Δz_surface = Oceananigans.Grids.zspacings(grid, Center(), Center(), Face(), Nz_grid)
    rate = restoring_rate / (Δz_surface * days)
    z_surface = Oceananigans.Grids.znode(Nz_grid, grid, Face()) - 1

    @inline surface_mask(x, y, z, t) = z >= z_surface
    FS = DatasetRestoring(Smetadata, arch; rate, mask=surface_mask, time_indices_in_memory=12)

    # Re-create ocean with forcing
    ocean = half_degree_tripolar_ocean(arch;
                                       Nz, depth,
                                       κ_skew, κ_symmetric,
                                       tracer_advection = WENO(order=7, minimum_buffer_upwind_order=3),
                                       forcing = (; S = FS))

    # Initial conditions from EN4
    set!(ocean.model, T=Metadatum(:temperature; dir=restoring_dir, dataset=EN4Monthly(), date=start_date),
                      S=Metadatum(:salinity;    dir=restoring_dir, dataset=EN4Monthly(), date=start_date))

    return ocean
end

# ----- ORCA dispatch -----

function _build_ocean(::Val{:orca}, arch;
                      Nz, depth, κ_skew, κ_symmetric,
                      restoring_dir, restoring_rate, start_date, end_date)

    ocean = orca_ocean(arch; Nz, depth, κ_skew, κ_symmetric)

    grid = ocean.model.grid

    # Surface salinity restoring to WOA
    Smetadata = Metadata(:salinity; dir=restoring_dir, dataset=WOAMonthly(),
                         start_date, end_date)

    Nz_grid = size(grid, 3)
    Δz_surface = Oceananigans.Grids.zspacings(grid, Center(), Center(), Face(), Nz_grid)
    rate = restoring_rate / (Δz_surface * days)
    z_surface = Oceananigans.Grids.znode(Nz_grid, grid, Face()) - 1

    @inline surface_mask(x, y, z, t) = z >= z_surface
    FS = DatasetRestoring(Smetadata, arch; rate, mask=surface_mask, time_indices_in_memory=12)

    ocean = orca_ocean(arch; Nz, depth, κ_skew, κ_symmetric,
                       forcing = (; S = FS))

    # Initial conditions from WOA
    set!(ocean.model, T=Metadatum(:temperature; dir=restoring_dir, dataset=WOAMonthly(), date=start_date),
                      S=Metadatum(:salinity;    dir=restoring_dir, dataset=WOAMonthly(), date=start_date))

    return ocean
end

# ----- Sea ice dispatch -----

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

# ----- Restart -----

function _load_restart!(coupled, atmosphere, ocean, sea_ice, checkpoint_path)
    # Load state from JLD2 checkpoint
    file = jldopen(checkpoint_path)

    # Ocean state
    set!(ocean.model, u = file["uo"],
                      v = file["vo"],
                      T = file["To"],
                      S = file["So"])

    if haskey(file, "eo")
        set!(ocean.model, e = file["eo"])
    end

    interior(ocean.model.free_surface.displacement) .= file["ηo"]
    interior(ocean.model.velocities.w) .= file["wo"]

    # Sea ice state
    interior(sea_ice.model.ice_thickness) .= file["hi"]
    interior(sea_ice.model.ice_concentration) .= file["ℵi"]
    interior(sea_ice.model.velocities.u) .= file["ui"]
    interior(sea_ice.model.velocities.v) .= file["vi"]

    # Sync clocks
    clock = file["clock"]
    ocean.model.clock.time = clock.time
    ocean.model.clock.iteration = clock.iteration
    sea_ice.model.clock.time = clock.time
    sea_ice.model.clock.iteration = clock.iteration
    coupled.clock.time = clock.time
    coupled.clock.iteration = clock.iteration

    # Sync atmosphere
    atmosphere.clock.time = clock.time
    atmosphere.clock.iteration = clock.iteration
    time_step!(atmosphere, 0)
    atmosphere.clock.iteration -= 1

    close(file)

    return nothing
end
```

**Important implementation notes:**
- The `_build_ocean` functions create the ocean twice — first to get the grid for computing restoring parameters, then with forcing. This mirrors the current `omip_defaults.jl` approach. At implementation time, check whether there's a way to add forcing after construction to avoid the double build. If `ocean_simulation` supports a `forcing` kwarg, building once should be sufficient — in that case, construct the grid first, compute restoring, then build the ocean once.
- The `@inline surface_mask` function captures `z_surface` from the enclosing scope. Verify this works with GPU compilation; if not, use a `parameters` approach.
- The exact `Metadata` vs `Metadatum` API and their keyword arguments should be verified against the current NumericalEarth API at implementation time.

### Step 3: Commit

```
feat(OMIPConfigurations): implement omip_simulation for half_degree and orca
```

---

## Task 4: Implement `add_omip_diagnostics!`

**Files:**
- Write: `src/OMIPConfigurations/diagnostics.jl`

This is the core of the diagnostics pipeline. It creates 4 output writers (surface, fields, scalars, checkpoint) and attaches them to the simulation.

### Step 1: Implement `diagnostics.jl`

```julia
using Oceananigans.OutputWriters: NetCDFOutputWriter, Checkpointer, JLD2OutputWriter
using Oceananigans.AbstractOperations: KernelFunctionOperation

"""
    add_omip_diagnostics!(simulation;
                          surface_averaging_interval = 15days,
                          field_averaging_interval = 15days,
                          checkpoint_interval = 90days,
                          output_dir = ".",
                          filename_prefix = "omip",
                          transports = false)

Attach standard OMIP diagnostics to a coupled ocean-sea-ice `simulation`.

Creates the following output writers:
- `surface`: 2D surface fields and fluxes (NetCDF), averaged over `surface_averaging_interval`
- `fields`: 3D prognostic fields (NetCDF), averaged over `field_averaging_interval`
- `scalars`: zonal/global mean quantities (NetCDF), averaged over `field_averaging_interval`
- `checkpointer`: full model state (JLD2), saved every `checkpoint_interval`

If `transports=true`, also adds MOC, heat transports, and strait transports.
"""
function add_omip_diagnostics!(simulation;
                               surface_averaging_interval = 15days,
                               field_averaging_interval = 15days,
                               checkpoint_interval = 90days,
                               output_dir = ".",
                               filename_prefix = "omip",
                               transports = false)

    model = simulation.model
    ocean = model.ocean
    sea_ice = model.sea_ice
    grid = ocean.model.grid

    # ============================
    # 1. Surface diagnostics (2D)
    # ============================

    Nz = size(grid, 3)

    # --- Ocean surface fields ---
    T, S = ocean.model.tracers.T, ocean.model.tracers.S
    u, v, w = ocean.model.velocities
    η = ocean.model.free_surface.displacement

    tos  = Field(T;    indices=(:, :, Nz))   # SST
    sos  = Field(S;    indices=(:, :, Nz))   # SSS
    zos  = Field(η)                           # SSH (already 2D)
    uo_s = Field(u;    indices=(:, :, Nz))   # surface u
    vo_s = Field(v;    indices=(:, :, Nz))   # surface v

    # Squared fields for variance computation
    tossq = Field(T^2; indices=(:, :, Nz))
    sossq = Field(S^2; indices=(:, :, Nz))
    zossq = Field(η^2)

    # Mixed layer depth
    mld = MixedLayerDepthField(ocean.model.buoyancy, grid, ocean.model.tracers)

    # --- Fluxes ---
    τx = model.interfaces.net_fluxes.ocean.u         # wind stress x
    τy = model.interfaces.net_fluxes.ocean.v         # wind stress y
    JT = model.interfaces.net_fluxes.ocean.T         # net heat flux (hfds)
    Js = model.interfaces.net_fluxes.ocean.S         # net salt/freshwater flux (wfo)
    Qc = model.interfaces.atmosphere_ocean_interface.fluxes.sensible_heat  # hfss
    Qv = model.interfaces.atmosphere_ocean_interface.fluxes.latent_heat    # hfls

    # --- Sea ice fields ---
    hi = sea_ice.model.ice_thickness
    ℵi = sea_ice.model.ice_concentration
    ui, vi = sea_ice.model.velocities

    surface_outputs = Dict(
        "tos"   => tos,
        "sos"   => sos,
        "zos"   => zos,
        "uo_s"  => uo_s,
        "vo_s"  => vo_s,
        "tossq" => tossq,
        "sossq" => sossq,
        "zossq" => zossq,
        "mlotst" => mld,
        "tauuo" => τx,
        "tauvo" => τy,
        "hfds"  => JT,
        "wfo"   => Js,
        "hfss"  => Qc,
        "hfls"  => Qv,
        "siconc"   => ℵi,
        "sithick"  => hi,
        "siu"      => ui,
        "siv"      => vi,
    )

    # Add sea ice surface temperature if available
    if hasproperty(sea_ice.model, :ice_thermodynamics) &&
       hasproperty(sea_ice.model.ice_thermodynamics, :top_surface_temperature)
        surface_outputs["sitemptop"] = sea_ice.model.ice_thermodynamics.top_surface_temperature
    end

    simulation.output_writers[:surface] = NetCDFOutputWriter(ocean.model, surface_outputs;
        filename = joinpath(output_dir, filename_prefix * "_surface.nc"),
        schedule = AveragedTimeInterval(surface_averaging_interval),
        overwrite_existing = true)

    # ============================
    # 2. 3D field diagnostics
    # ============================

    field_outputs = Dict(
        "thetao" => T,
        "so"     => S,
        "uo"     => u,
        "vo"     => v,
        "wo"     => w,
    )

    # TKE from CATKE (if present)
    if haskey(ocean.model.tracers, :e)
        field_outputs["tke"] = ocean.model.tracers.e
    end

    simulation.output_writers[:fields] = NetCDFOutputWriter(ocean.model, field_outputs;
        filename = joinpath(output_dir, filename_prefix * "_fields.nc"),
        schedule = AveragedTimeInterval(field_averaging_interval),
        overwrite_existing = true)

    # ============================
    # 3. Scalar / zonal means
    # ============================
    # NOTE: Zonal means require online reduction infrastructure.
    # Implementation depends on Oceananigans' `Average` and `Reduction` support.
    # The exact approach should be verified at implementation time.
    # Options:
    #   a) Use Oceananigans' built-in `Average(field, dims=1)` if it works
    #      with immersed grids and NetCDF output
    #   b) Use a callback that computes zonal means manually (like
    #      `compute_zonal_averages` in Diagnostics/report_computations.jl)
    #      and writes to a separate file
    #   c) Use `Reduction` operations
    #
    # For global scalars, use `Reduction` or a callback:

    # Global mean scalars via callback
    # This writes a small NetCDF file with global/zonal means at each output time.
    # Implementation skeleton:

    # scalar_callback = OMIPScalarDiagnostics(grid, T, S, ocean;
    #     filename = joinpath(output_dir, filename_prefix * "_scalars.nc"),
    #     schedule = AveragedTimeInterval(field_averaging_interval))
    # add_callback!(simulation, scalar_callback, TimeInterval(field_averaging_interval))

    # TODO: Implement OMIPScalarDiagnostics struct + callback.
    # This needs:
    #   - Global volume-weighted mean T, S (thetaoga, soga)
    #   - Global area-weighted mean SST, SSS (tosga, sosga)
    #   - Total ocean mass and volume (masso, volo)
    #   - Zonal mean T(y,z) and S(y,z)
    # Use the existing `compute_zonal_averages` from Diagnostics module as reference.

    # ============================
    # 4. Checkpointer
    # ============================

    simulation.output_writers[:checkpointer] = Checkpointer(model;
        schedule = TimeInterval(checkpoint_interval),
        prefix = joinpath(output_dir, filename_prefix * "_checkpoint"),
        cleanup = false)

    # ============================
    # 5. Progress callback
    # ============================

    progress = Progress()
    add_callback!(simulation, progress, IterationInterval(10))

    # ============================
    # 6. Transport diagnostics (opt-in)
    # ============================

    if transports
        add_omip_transport_diagnostics!(simulation;
            output_dir, filename_prefix, field_averaging_interval)
    end

    return nothing
end
```

**Implementation notes:**
- The exact NetCDF output writer API depends on Oceananigans version. At implementation time, verify:
  - Whether `NetCDFOutputWriter` accepts a `Dict` or needs a `NamedTuple`
  - Whether it can write `Field` slices (indices) and `AbstractOperations` (T^2)
  - Whether `AveragedTimeInterval` is the correct schedule type
- The flux field access pattern (`model.interfaces.net_fluxes...`) is taken from `omip_defaults.jl:173-181`. Verify these paths still work with the current NumericalEarth coupled model API.
- The scalar/zonal mean diagnostics are the most complex piece — they may need a custom callback rather than an output writer. The `compute_zonal_averages` function in `Diagnostics/report_computations.jl:67-122` provides the algorithm but computes on CPU snapshots. For online averaging, either use Oceananigans `Average` reductions or implement a custom periodic callback.
- Sea ice fields may need wrapping in `Field(...)` if they aren't already `AbstractField` objects compatible with the output writer.

### Step 2: Commit

```
feat(OMIPConfigurations): implement add_omip_diagnostics! with surface, 3D, and checkpoint output
```

---

## Task 5: Implement transport diagnostics (opt-in)

**Files:**
- Write: `src/OMIPConfigurations/transport_diagnostics.jl`

### Step 1: Implement `transport_diagnostics.jl`

```julia
"""
    add_omip_transport_diagnostics!(simulation; output_dir, filename_prefix, field_averaging_interval)

Add opt-in transport diagnostics: MOC streamfunction, meridional heat/salt
transports, and key strait mass transports.

These are Priority 1 in the OMIP protocol but computationally expensive.
"""
function add_omip_transport_diagnostics!(simulation;
                                         output_dir = ".",
                                         filename_prefix = "omip",
                                         field_averaging_interval = 15days)

    # TODO: Implementation depends on available infrastructure in
    # ClimaOcean.Diagnostics (MeridionalStreamfunction, compute_amoc, etc.)
    # and basin mask definitions.
    #
    # Required components:
    #
    # 1. MOC streamfunction (msftmyz):
    #    - Requires basin masks (Atlantic-Arctic, Indo-Pacific, Global)
    #    - Integrate v*dz across longitude for each basin
    #    - Output: latitude × depth × basin
    #    - Use `MeridionalStreamfunction` from Diagnostics module if available
    #
    # 2. Meridional heat transport (hfbasin):
    #    - ∫∫ ρ₀ cₚ v T dx dz for each basin
    #    - Output: latitude × basin
    #
    # 3. Meridional salt transport:
    #    - ∫∫ ρ₀ v S dx dz for each basin
    #    - Output: latitude × basin
    #
    # 4. Strait mass transports (mfo):
    #    - Drake Passage, Bering Strait, Indonesian Throughflow,
    #      Florida-Bahamas, Mozambique Channel, etc.
    #    - Requires section definitions (pairs of lat/lon endpoints)
    #    - Integrate ρ₀ u·n dA across each section
    #
    # Implementation approach:
    #    - Use a TimeInterval callback that computes these at each output time
    #    - Write results to a separate NetCDF file
    #    - Basin masks can be computed from bathymetry + flood fill
    #    - Section transports via interpolation to section coordinates
    #
    # This is a separate development effort. For now, log a warning:

    @warn "Transport diagnostics (transports=true) not yet implemented. " *
          "MOC, heat transports, and strait transports will be added in a future release."

    return nothing
end
```

### Step 2: Commit

```
feat(OMIPConfigurations): add transport diagnostics stub (opt-in, not yet implemented)
```

---

## Task 6: Implement scalar/zonal mean diagnostics

**Files:**
- Modify: `src/OMIPConfigurations/diagnostics.jl`

This task fills in the scalar diagnostics that were left as TODO in Task 4.

### Step 1: Implement `OMIPScalarDiagnostics` callback

Add to `diagnostics.jl`:

```julia
using NCDatasets

"""
    OMIPScalarCallback

Callback that computes and appends zonal means and global scalars to a NetCDF file
at each invocation.
"""
mutable struct OMIPScalarCallback{G, T, S, F}
    grid :: G
    T_field :: T
    S_field :: S
    filepath :: F
    initialized :: Bool
end

function (cb::OMIPScalarCallback)(simulation)
    grid = cb.grid
    T = cb.T_field
    S = cb.S_field
    t = time(simulation)

    # Compute zonal means using existing infrastructure
    T̄, S̄, φ, z = compute_zonal_averages(grid, T, S)

    # Compute global means
    Nz = size(grid, 3)
    SST = Array(interior(T, :, :, Nz))
    SSS = Array(interior(S, :, :, Nz))

    tosga = mean(filter(!isnan, SST))
    sosga = mean(filter(!isnan, SSS))
    thetaoga = mean(filter(!isnan, Array(interior(T))))
    soga = mean(filter(!isnan, Array(interior(S))))

    # Write to NetCDF
    if !cb.initialized
        _initialize_scalar_nc(cb.filepath, φ, z)
        cb.initialized = true
    end
    _append_scalar_nc(cb.filepath, t, T̄, S̄, tosga, sosga, thetaoga, soga)

    return nothing
end

function _initialize_scalar_nc(filepath, φ, z)
    ds = NCDataset(filepath, "c")

    defDim(ds, "latitude", length(φ))
    defDim(ds, "depth", length(z))
    defDim(ds, "time", Inf)  # unlimited

    defVar(ds, "latitude", Float64, ("latitude",))[:] = φ
    defVar(ds, "depth", Float64, ("depth",))[:] = z
    defVar(ds, "time", Float64, ("time",))

    defVar(ds, "T_zonal_mean", Float32, ("latitude", "depth", "time"))
    defVar(ds, "S_zonal_mean", Float32, ("latitude", "depth", "time"))
    defVar(ds, "tosga", Float32, ("time",))
    defVar(ds, "sosga", Float32, ("time",))
    defVar(ds, "thetaoga", Float32, ("time",))
    defVar(ds, "soga", Float32, ("time",))

    close(ds)
end

function _append_scalar_nc(filepath, t, T̄, S̄, tosga, sosga, thetaoga, soga)
    ds = NCDataset(filepath, "a")

    n = length(ds["time"]) + 1
    ds["time"][n] = t
    ds["T_zonal_mean"][:, :, n] = Float32.(T̄)
    ds["S_zonal_mean"][:, :, n] = Float32.(S̄)
    ds["tosga"][n] = Float32(tosga)
    ds["sosga"][n] = Float32(sosga)
    ds["thetaoga"][n] = Float32(thetaoga)
    ds["soga"][n] = Float32(soga)

    close(ds)
end
```

Then in `add_omip_diagnostics!`, replace the TODO block with:

```julia
    # Scalar / zonal mean diagnostics
    scalar_cb = OMIPScalarCallback(grid, T, S,
        joinpath(output_dir, filename_prefix * "_scalars.nc"),
        false)

    add_callback!(simulation, scalar_cb, TimeInterval(field_averaging_interval))
```

**Implementation notes:**
- `compute_zonal_averages` is imported from `Diagnostics.report_computations` — it handles both LatLon and Tripolar grids
- The global means computed here are **snapshot** means, not time-averaged means. For proper time-averaged global means, accumulation logic would be needed. For a first implementation, snapshot means at the callback time are acceptable — the 3D fields in the `fields` output writer are already time-averaged, so users can compute global means from those if needed.
- NCDatasets.jl may already be a dependency via NumericalEarth — check `Project.toml`. If not, it needs to be added.

### Step 2: Commit

```
feat(OMIPConfigurations): implement zonal mean and global scalar diagnostics
```

---

## Task 7: Update experiment scripts in NumericalEarth.jl

**Files:**
- Modify: `NumericalEarth.jl/experiments/omips/half_degree_omip/half_degree_omip.jl`
- Modify: `NumericalEarth.jl/experiments/omips/orca_omip/orca_omip.jl`

### Step 1: Simplify `half_degree_omip.jl`

Replace the entire file with:

```julia
using ClimaOcean.OMIPConfigurations

sim = omip_simulation(:half_degree;
    restart = "halfdegree_iteration18432_checkpoint.jld2",
    forcing_dir = "forcing_data",
    restoring_dir = "climatology",
    filename_prefix = "halfdegree",
    Δt = 30minutes,
    stop_time = Inf,
)

run!(sim)
```

### Step 2: Simplify `orca_omip.jl`

Replace with:

```julia
using ClimaOcean.OMIPConfigurations

sim = omip_simulation(:orca;
    arch = GPU(),
    restoring_dir = "climatology",
    filename_prefix = "orca",
    Δt = 30minutes,
    stop_time = 300 * 365days,
)

run!(sim)
```

### Step 3: Commit

```
refactor(experiments): simplify OMIP scripts to use ClimaOcean.OMIPConfigurations
```

---

## Task 8: Testing

**Files:**
- Create: `test/test_omip_configurations.jl`
- Modify: `test/runtests.jl` (add include)

### Step 1: Write unit tests

```julia
using Test
using ClimaOcean
using ClimaOcean.OMIPConfigurations
using Oceananigans
using Oceananigans.Units

@testset "OceanConfigurations parameter API" begin
    # Test that physical parameters propagate correctly
    # Use a small grid to keep tests fast
    # NOTE: These tests require data downloads; mark as integration tests
    #       or mock the data layer

    @testset "half_degree_tripolar_ocean accepts Nz kwarg" begin
        # Verify the function signature accepts Nz without error
        # (actual grid construction requires bathymetry data)
        m = methods(half_degree_tripolar_ocean)
        @test length(m) > 0
    end

    @testset "default_half_degree_closure accepts kwargs" begin
        closure = ClimaOcean.OceanConfigurations.default_half_degree_closure(;
            κ_skew = 300,
            κ_symmetric = 100,
            biharmonic_timescale = 20days,
        )
        @test length(closure) == 4
        # Check eddy closure has correct parameters
        @test closure[2].κ_skew == 300
        @test closure[2].κ_symmetric == 100
    end
end

@testset "add_omip_diagnostics! API" begin
    # Test that the function exists and has the expected signature
    m = methods(add_omip_diagnostics!)
    @test length(m) > 0
end
```

**Note:** Full integration tests (actually building simulations) require bathymetry/forcing data downloads and are slow. These should be marked as integration tests or run in CI with data caching. The unit tests above verify API signatures and parameter propagation without requiring data.

### Step 2: Commit

```
test: add basic tests for OMIPConfigurations API
```

---

## Summary of files changed/created

### New files (in `ClimaOcean.jl`):
1. `src/OMIPConfigurations/OMIPConfigurations.jl` — module definition
2. `src/OMIPConfigurations/omip_simulation.jl` — main entry point + config dispatches
3. `src/OMIPConfigurations/atmosphere.jl` — shared JRA55 atmosphere setup
4. `src/OMIPConfigurations/diagnostics.jl` — `add_omip_diagnostics!` + scalar callback
5. `src/OMIPConfigurations/transport_diagnostics.jl` — opt-in transport diagnostics (stub)
6. `test/test_omip_configurations.jl` — tests

### Modified files (in `ClimaOcean.jl`):
1. `src/ClimaOcean.jl` — add module include + exports
2. `src/OceanConfigurations/OceanConfigurations.jl` — update `vertical_coordinate` signature
3. `src/OceanConfigurations/half_degree_tripolar.jl` — parameter-based API
4. `src/OceanConfigurations/orca.jl` — parameter-based API
5. `src/OceanConfigurations/one_degree_tripolar.jl` — parameter-based API (consistency)

### Modified files (in `NumericalEarth.jl`):
1. `experiments/omips/half_degree_omip/half_degree_omip.jl` — simplified to 10 lines
2. `experiments/omips/orca_omip/orca_omip.jl` — simplified to 10 lines

### Files that become obsolete:
1. `NumericalEarth.jl/experiments/omips/omip_defaults.jl` — logic migrated to OMIPConfigurations
2. `NumericalEarth.jl/experiments/omips/initialize_omip.jl` — restart logic migrated to `_load_restart!`

---

## Execution order

Tasks 1-6 are sequential (each builds on the previous). Task 7 (experiment scripts) depends on Tasks 1-4. Task 8 (testing) can be done alongside each task.

Recommended order: **1 → 2 → 3 → 4 → 6 → 5 → 7 → 8**

(Transport diagnostics (Task 5) is a stub, so it can be done anytime. Scalar diagnostics (Task 6) should come right after the main diagnostics (Task 4) since it fills in the TODO.)
