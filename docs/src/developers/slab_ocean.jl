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
# To integrate a slab ocean into `OceanSeaIceModel`, we need to extend several methods from ClimaOcean's `OceanSeaIceModels.jl` and
# `InterfaceComputations.jl` modules as well as Oceananigans' TimeSteppers.jl module
#
# 1. **ComponentExchanger**: Defines how to extract state from the component for flux computations
# 2. **interpolate_state!**: Interpolates component state onto the exchange grid
# 3. **update_net_fluxes!**: Computes net fluxes and applies them to the component
# 4. **reference_density**: of the ocean component
# 5  **heat_capacity**: of the ocean component
# 6. **ocean_salinity**: returns the entire salinity field of the ocean component
# 7. **ocean_temperature**: returns the entire temperature field of the ocean component
# 8. **ocean_surface_salinity**: returns the salinity at the surface of the ocean component
# 8. **net_fluxes**: Returns the flux fields that will be updated by the coupling system
# 6. **time_step!**: Advances the component forward in time
#
# ## Step 1: Define the Slab Ocean Type
#
# First, we define a struct to hold the slab ocean state and its constructor. We also define a summary and show method for the slab ocean type.

using Oceananigans, Base

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

Base.summary(slab_ocean::SlabOcean) = "SlabOcean with Depth=$(slab_ocean.grid.Lz)"
Base.show(io::IO, slab_ocean::SlabOcean) = print(io, Base.summary(slab_ocean))

# ## Step 2: Extend the InterfaceComputations.jl module
#
# The `ComponentExchanger` type contains the state variables needed for flux computations.
# These are the ocean surface state on the `exchange_grid` and the "regridder" to interpolate data from the ocean onto the `exchange_grid`.
# Here, we assume that the ocean is on the same grid as the exchange grid, so that the regridder is ``nothing'' and the state variables are the same as the ocean surface state.
# The flux computation requires also ocean surface velocities to compute the turbulent fluxes. In this case, we use dummy zero velocity fields since the slab ocean has no dynamics.

using ClimaOcean.OceanSeaIceModels.InterfaceComputations

function InterfaceComputations.ComponentExchanger(slab_ocean::SlabOcean, exchange_grid)
    T = slab_ocean.temperature
    S = slab_ocean.salinity
    u = Oceananigans.Fields.ZeroField()
    v = Oceananigans.Fields.ZeroField()
    return InterfaceComputations.ComponentExchanger((; u, v, T, S), nothing)
end

# The `net_fluxes` function returns the flux fields that will be updated by the coupling system. 
# For a slab ocean, we need to return the containers for the temperature and freshwater (salinity) fluxes
# stored in the `slab_ocean.fluxes` field as well as the dummy stress fields which will be unused since the slab ocean has no dynamics.

function InterfaceComputations.net_fluxes(slab_ocean::SlabOcean)
    grid = slab_ocean.grid
    τx = Field{Center, Center, Nothing}(grid)
    τy = Field{Center, Center, Nothing}(grid)
    return merge(slab_ocean.fluxes, (; u=τx, v=τy))
end 

# ## Step 3: Extend the OceanSeaIceModels.jl module
#
# In the OceanSeaIceModels.jl module, we define the thermodynamic properties of the ocean component as well as all the helper functions 
# needed to retrieve the ocean state and surface state.

using ClimaOcean.OceanSeaIceModels

OceanSeaIceModels.reference_density(::SlabOcean) = 1025.0
OceanSeaIceModels.heat_capacity(::SlabOcean) = 3990.0
OceanSeaIceModels.ocean_surface_salinity(slab_ocean::SlabOcean) = slab_ocean.salinity
OceanSeaIceModels.ocean_salinity(slab_ocean::SlabOcean) = slab_ocean.salinity
OceanSeaIceModels.ocean_temperature(slab_ocean::SlabOcean) = slab_ocean.temperature

# The `update_net_fluxes!` function computes net fluxes and applies them to previously defined ``net_fluxes'' containers.
# These will be used to update the ocean state in the `time_step!` method. In this case, we can use the `update_net_ocean_fluxes!` function from the Oceans.jl module.

using ClimaOcean.Oceans

OceanSeaIceModels.update_net_fluxes!(coupled_model, slab_ocean::SlabOcean) = 
    Oceans.update_net_ocean_fluxes!(coupled_model, slab_ocean, slab_ocean.grid)

# ## Step 4: Extend the TimeSteppers.jl module
#
# The `time_step!` method is called by the coupled model to advance the component forward in time. 
# For a slab ocean, this method advances temperature and salinity through the computed fluxes:

import Oceananigans.TimeSteppers: time_step!
using Oceananigans.TimeSteppers: tick!

function time_step!(slab_ocean::SlabOcean, Δt)
    tick!(slab_ocean.clock, Δt)    
    parent(slab_ocean.temperature) .+= parent(slab_ocean.fluxes.T) .* Δt ./ size(slab_ocean.grid, 3)
    parent(slab_ocean.salinity)    .+= parent(slab_ocean.fluxes.S) .* Δt ./ size(slab_ocean.grid, 3)
    return nothing
end

# ## Step 5: Complete Example: Coupling Slab Ocean with JRA55 and Sea Ice
#
# Here's a complete example showing how to use the slab ocean in a coupled simulation:

using ClimaOcean
using Oceananigans
using Oceananigans.Units
using Dates

arch = CPU()
grid = LatitudeLongitudeGrid(arch,
                             size = (720, 360, 1),  # 1 degree resolution, single vertical level
                             longitude = (0, 360),
                             latitude = (-90, 90),
                             z = (-50, 0))  # Single layer from -50m to 0m (surface)

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

slab_ocean = SlabOcean(grid)
set!(slab_ocean.temperature, Metadatum(:temperature, dataset=ECCO4Monthly()))
set!(slab_ocean.salinity,    Metadatum(:salinity,    dataset=ECCO4Monthly()))

atmosphere = ClimaOcean.JRA55PrescribedAtmosphere(arch)

sea_ice = ClimaOcean.sea_ice_simulation(grid, slab_ocean)
set!(sea_ice.model, h=1, ℵ=1)

interfaces = ComponentInterfaces(atmosphere, slab_ocean, sea_ice; exchange_grid=grid)
coupled_model = ClimaOcean.OceanSeaIceModel(slab_ocean, sea_ice; atmosphere, interfaces)

simulation = Simulation(coupled_model, Δt = 30minutes, stop_time = 100days)
run!(simulation)

using CairoMakie

fig = Figure()
axT = Axis(fig[1, 1])
axS = Axis(fig[1, 2])
axh = Axis(fig[2, 1])
axℵ = Axis(fig[2, 2])
heatmap!(axT, slab_ocean.temperature)
heatmap!(axS, slab_ocean.salinity)
heatmap!(axh, sea_ice.model.ice_thickness)
heatmap!(axℵ, sea_ice.model.ice_concentration)
Colorbar(fig[2, 1], hmT, vertical=false)
Colorbar(fig[2, 2], hmS, vertical=false)

display(fig)

# ## Summary
#
# The slab ocean example demonstrates how a simplified component can be integrated into the coupling framework while maintaining compatibility 
# with the existing atmosphere and sea ice components. 
# The key insight is that the coupling system handles flux computations generically; 
# your component just needs to provide the right interface methods.
