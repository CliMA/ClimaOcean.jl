# # Implementing a Slab Ocean Component
#
# This tutorial demonstrates how to implement a new ocean component for `OceanSeaIceModel`: a **slab ocean model**. A slab ocean represents the ocean as a single well-mixed layer with a fixed depth, making it computationally efficient while still capturing the essential thermodynamic coupling between the ocean, atmosphere, and sea ice. This approach is commonly used in climate sensitivity studies and is described in detail in [Garuba et al. (2024)](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2023MS003910).
#
# ## Overview
#
# A slab ocean model simplifies the full 3D ocean by:
# - Representing the ocean as a single temperature (and optionally salinity) per horizontal grid point
# - Using a fixed "slab depth" that determines the heat capacity
# - Evolving temperature based on net heat fluxes from the atmosphere and sea ice
# - Having no ocean dynamics (no velocities)
#
# To integrate a slab ocean into `OceanSeaIceModel`, we need to extend several methods from `OceanSeaIceModels.jl` and `InterfaceComputations.jl`:
#
# 1. **ComponentExchanger**: Defines how to extract state from the component for flux computations
# 2. **interpolate_state!**: Interpolates component state onto the exchange grid
# 3. **update_net_fluxes!**: Computes net fluxes and applies them to the component
# 4. **reference_density** and **heat_capacity**: Provide thermodynamic properties
# 5. **net_fluxes**: Returns the flux fields that will be updated by the coupling system
#
# ## Step 1: Define the Slab Ocean Type
#
# First, we define a struct to hold the slab ocean state:

using Oceananigans

struct SlabOcean{G, C, T, S, F}
    grid :: G
    clock :: C
    temperature :: T   
    salinity :: S
    fluxes :: F
end

function SlabOcean(grid; clock =  Clock(time=0))
    temperature = CenterField(grid)
    salinity = CenterField(grid)
    fluxes = (T = CenterField(grid),
              S = CenterField(grid))

    return SlabOcean(grid, clock, temperature, salinity, fluxes)
end

# ## Step 2: Implement ComponentExchanger
#
# The `ComponentExchanger` extracts the state variables needed for flux computations.
# For a slab ocean, we need to provide temperature and salinity (and dummy velocity fields since the slab has no dynamics):

using ClimaOcean.OceanSeaIceModels.InterfaceComputations

function InterfaceComputations.ComponentExchanger(slab_ocean::SlabOcean, grid)
    T = slab_ocean.temperature
    S = slab_ocean.salinity
    u = Oceananigans.Fields.ZeroField()
    v = Oceananigans.Fields.ZeroField()
    return ComponentExchanger((; u, v, T, S), nothing)
end

# ## Step 4: Implement reference_density and heat_capacity
#
# These functions provide thermodynamic properties needed for flux computations as well as...

using ClimaOcean.OceanSeaIceModels

OceanSeaIceModels.reference_density(::SlabOcean) = 1025.0
OceanSeaIceModels.heat_capacity(::SlabOcean) = 3990.0
OceanSeaIceModels.ocean_surface_salinity(slab_ocean::SlabOcean) = slab_ocean.salinity
OceanSeaIceModels.ocean_salinity(slab_ocean::SlabOcean) = slab_ocean.salinity
OceanSeaIceModels.ocean_temperature(slab_ocean::SlabOcean) = slab_ocean.temperature

# ## Step 5: Implement net_fluxes
#
# The `net_fluxes` function returns the flux fields that will be updated by the coupling system. For a slab ocean, we need heat and freshwater (salinity) fluxes:

function InterfaceComputations.net_fluxes(slab_ocean::SlabOcean)
    grid = slab_ocean.grid
    τx = Field{Center, Center, Nothing}(grid)
    τy = Field{Center, Center, Nothing}(grid)
    return merge(slab_ocean.fluxes, (; u=τx, v=τy))
end

# ## Step 6: Implement update_net_fluxes!
#
# This is the core function that computes net fluxes and applies them to update the slab ocean temperature. Unlike a full ocean model where fluxes are applied as boundary conditions and the time-stepper handles updates, a slab ocean needs to update its temperature directly within `update_net_fluxes!`. The flux assembly follows the same pattern as the full ocean model:

using ClimaOcean.Oceans

OceanSeaIceModels.update_net_fluxes!(coupled_model, slab_ocean::SlabOcean) = 
    Oceans.update_net_ocean_fluxes!(coupled_model, slab_ocean, slab_ocean.grid)

# ## Step 7: Implement time_step!
#
# The `time_step!` method is called by the coupled model to advance the component forward in time. For a slab ocean, this method advances temperature and salinity through the computed fluxes:

import Oceananigans.TimeSteppers: time_step!
using Oceananigans.TimeSteppers: tick!

function time_step!(slab_ocean::SlabOcean, Δt)
    tick!(slab_ocean.clock, Δt)    
    parent(slab_ocean.temperature) .+= parent(slab_ocean.fluxes.T) .* Δt ./ size(slab_ocean.grid, 3)
    parent(slab_ocean.salinity)    .+= parent(slab_ocean.fluxes.S) .* Δt ./ size(slab_ocean.grid, 3)
    return nothing
end

# ## Step 7: Complete Example: Coupling Slab Ocean with JRA55 and Sea Ice
#
# Here's a complete example showing how to use the slab ocean in a coupled simulation:

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
sea_ice = ClimaOcean.sea_ice_simulation(grid, slab_ocean)
set!(sea_ice.model, h=1, ℵ=1)

# Create coupled model
interfaces = ComponentInterfaces(atmosphere, ocean, sea_ice; exchange_grid=grid)
coupled_model = ClimaOcean.OceanSeaIceModel(slab_ocean, sea_ice; atmosphere, interfaces)

# Create and run simulation
simulation = Simulation(coupled_model, Δt = 30minutes, stop_time = 100days)
run!(simulation)

# ## Summary
#
# To implement a new ocean component for `OceanSeaIceModel`, you need to extend:
#
# 1. **`ComponentExchanger(component, grid)`**: Extract state variables for flux computations
# 2. **`interpolate_state!(exchanger, grid, component, coupled_model)`**: Interpolate component state to exchange grid
# 3. **`update_net_fluxes!(coupled_model, component)`**: Compute and apply net fluxes
# 4. **`reference_density(component)`**: Return reference density (kg m⁻³)
# 5. **`heat_capacity(component)`**: Return heat capacity (J m⁻² K⁻¹ for 2D, J m⁻³ K⁻¹ for 3D)
# 6. **`net_fluxes(component)`**: Return NamedTuple of flux fields `(; u, v, T, S, ...)`
#
# The slab ocean example demonstrates how a simplified component can be integrated into the coupling framework while maintaining compatibility with the existing atmosphere and sea ice components. The key insight is that the coupling system handles flux computations generically; your component just needs to provide the right interface methods.
