# Implementing a Slab Ocean Component

This tutorial demonstrates how to implement a new ocean component for `OceanSeaIceModel`: a **slab ocean model**. A slab ocean represents the ocean as a single well-mixed layer with a fixed depth, making it computationally efficient while still capturing the essential thermodynamic coupling between the ocean, atmosphere, and sea ice. This approach is commonly used in climate sensitivity studies and is described in detail in [Garuba et al. (2024)](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2023MS003910).

## Overview

A slab ocean model simplifies the full 3D ocean by:
- Representing the ocean as a single temperature (and optionally salinity) per horizontal grid point
- Using a fixed "slab depth" that determines the heat capacity
- Evolving temperature based on net heat fluxes from the atmosphere and sea ice
- Having no ocean dynamics (no velocities)

To integrate a slab ocean into `OceanSeaIceModel`, we need to extend several methods from `OceanSeaIceModels.jl` and `InterfaceComputations.jl`:

1. **ComponentExchanger**: Defines how to extract state from the component for flux computations
2. **interpolate_state!**: Interpolates component state onto the exchange grid
3. **update_net_fluxes!**: Computes net fluxes and applies them to the component
4. **reference_density** and **heat_capacity**: Provide thermodynamic properties
5. **net_fluxes**: Returns the flux fields that will be updated by the coupling system

## Step 1: Define the Slab Ocean Type

First, we define a struct to hold the slab ocean state:

```julia
using Oceananigans

struct SlabOcean{G, C, T, S, F, FT}
    grid :: G
    clock :: C
    tracers :: T    # Field{Center, Center, Nothing} - single temperature per horizontal point
    salinity :: S   # Field{Center, Center, Nothing} - single salinity per horizontal point (optional)
    fluxes :: F
    depth :: FT     # Fixed depth of the slab (meters)
end

function SlabOcean(grid; slab_depth = 50)
    FT = eltype(grid)
    
    # Create 2D fields (no vertical dimension) for temperature and salinity
    tracers = (T = Field{Center, Center, Nothing}(grid),
               S = Field{Center, Center, Nothing}(grid))

    fluxes = deepcopy(tracers)

    return SlabOcean(grid, Clock(time=0), tracers, fluxes, convert(FT, slab_depth))
end
```

## Step 2: Implement ComponentExchanger

The `ComponentExchanger` extracts the state variables needed for flux computations.
For a slab ocean, we need to provide temperature and salinity (and dummy velocity fields since the slab has no dynamics):

```julia
import ClimaOcean.OceanSeaIceModels.InterfaceComputations: ComponentExchanger

function ComponentExchanger(slab_ocean::SlabOcean, grid)
    # The exchange grid should match the slab ocean grid horizontally
    # For a slab ocean, we create views/fields that represent the "surface" state
    
    # Temperature and salinity are already 2D (no vertical dimension)
    T = slab_ocean.tracers.T
    S = slab_ocean.tracers.S
    
    # Create zero velocity fields (slab ocean has no dynamics)
    # These are needed for compatibility with the flux computation interface
    u = ZeroField()
    v = ZeroField()

    # Return ComponentExchanger with state and no regridder (grids match)
    return ComponentExchanger((; u, v, T, S), nothing)
end
```

```julia
import ClimaOcean.SeaIceSimulations: ocean_surface_salinity

ocean_surface_salinity(slab_ocean::SlabOcean) = slab_ocean.tracers.T

```

## Step 4: Implement reference_density and heat_capacity

These functions provide thermodynamic properties needed for flux computations:

```julia
import ClimaOcean.OceanSeaIceModels: reference_density, heat_capacity

# Reference density for seawater (kg m⁻³)
reference_density(::SlabOcean) = 1025.0

# Specific heat capacity: 
heat_capacity(slab_ocean::SlabOcean) = 3990.0
```

## Step 5: Implement net_fluxes

The `net_fluxes` function returns the flux fields that will be updated by the coupling system. For a slab ocean, we need heat and freshwater (salinity) fluxes:

```julia
import ClimaOcean.OceanSeaIceModels.InterfaceComputations: net_fluxes

function net_fluxes(slab_ocean::SlabOcean)
    grid = slab_ocean.grid
    
    # Create flux fields for temperature (heat) and salinity (freshwater)
    # These will be filled by update_net_fluxes!
    Jᵀ = slab_ocean.fluxes.T # Temperature flux (K m s⁻¹)
    Jˢ = slab_ocean.fluxes.S # Salinity flux (psu m s⁻¹)
    
    # Dummy momentum fluxes for compatibilty with
    # the flux computation interface (slab ocean has no dynamics)
    τx = Field{Center, Center, Nothing}(grid)
    τy = Field{Center, Center, Nothing}(grid)
    
    return (; u=τx, v=τy, T=Jᵀ, S=Jˢ)
end
```

## Step 6: Implement update_net_fluxes!

This is the core function that computes net fluxes and applies them to update the slab ocean temperature. Unlike a full ocean model where fluxes are applied as boundary conditions and the time-stepper handles updates, a slab ocean needs to update its temperature directly within `update_net_fluxes!`. The flux assembly follows the same pattern as the full ocean model:

```julia
using ClimaOcean.Oceans
import ClimaOcean.OceanSeaIceModels: update_net_fluxes!

update_net_fluxes!(coupled_model, slab_ocean::SlabOcean) = Oceans.update_net_ocean_fluxes!(coupled_model, slab_ocean)
```

## Step 7: Implement time_step!

The `time_step!` method is called by the coupled model to advance the component forward in time. For a slab ocean, this method advances temperature and salinity through the computed fluxes:

```julia
import Oceananigans.TimeSteppers: time_step!
using Oceananigans.TimeSteppers: tick!

function time_step!(slab_ocean::SlabOcean, Δt; callbacks=[], compute_tendencies=true)
    # Tick the clock forward
    tick!(slab_ocean.clock, Δt)
    
    parent(slab_ocean.tracers.T) .+= parent(slab_ocean.Jᵀ) .* Δt ./ slab_ocean.depth
    parent(slab_ocean.tracers.S) .+= parent(slab_ocean.Jˢ) .* Δt ./ slab_ocean.depth

    return nothing
end
```

## Step 7: Complete Example: Coupling Slab Ocean with JRA55 and Sea Ice

Here's a complete example showing how to use the slab ocean in a coupled simulation:

```julia
using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Dates

# Create a horizontal grid with a single vertical level for slab ocean
arch = CPU()
grid = LatitudeLongitudeGrid(arch,
                             size = (720, 360, 1),  # 1 degree resolution, single vertical level
                             longitude = (0, 360),
                             latitude = (-90, 90),
                             z = (-50, 0))  # Single layer from -50m to 0m (surface)

# Create slab ocean
slab_ocean = SlabOcean(grid)
set!(slab_ocean.tracers.T, Metadatum(:temperature, dataset=ECCO4Monthly()))
set!(slab_ocean.tracers.S, Metadatum(:salinity,    dataset=ECCO4Monthly()))

# Create prescribed atmosphere
atmosphere = ClimaOcean.JRA55PrescribedAtmosphere(arch)

# Create sea ice simulation (on the same grid)
sea_ice = ClimaOcean.sea_ice_simulation(grid, ocean)
set!(sea_ice.model, h=1, ℵ=1)

# Create coupled model
coupled_model = ClimaOcean.OceanSeaIceModel(slab_ocean, sea_ice; atmosphere)

# Create and run simulation
simulation = Simulation(coupled_model, Δt = 30minutes, stop_time = 100days)
run!(simulation)
```

## Summary

To implement a new ocean component for `OceanSeaIceModel`, you need to extend:

1. **`ComponentExchanger(component, grid)`**: Extract state variables for flux computations
2. **`interpolate_state!(exchanger, grid, component, coupled_model)`**: Interpolate component state to exchange grid
3. **`update_net_fluxes!(coupled_model, component)`**: Compute and apply net fluxes
4. **`reference_density(component)`**: Return reference density (kg m⁻³)
5. **`heat_capacity(component)`**: Return heat capacity (J m⁻² K⁻¹ for 2D, J m⁻³ K⁻¹ for 3D)
6. **`net_fluxes(component)`**: Return NamedTuple of flux fields `(; u, v, T, S, ...)`

The slab ocean example demonstrates how a simplified component can be integrated into the coupling framework while maintaining compatibility with the existing atmosphere and sea ice components. The key insight is that the coupling system handles flux computations generically; your component just needs to provide the right interface methods.
