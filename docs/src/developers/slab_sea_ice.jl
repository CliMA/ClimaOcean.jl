# # Implementing a Slab Sea Ice Component
#
# This tutorial demonstrates how to implement a new sea ice component for `OceanSeaIceModel`: a **slab sea ice model**.
# A slab sea ice represents sea ice as thickness and concentration fields evolving based on heat fluxes,
# without the complexity of ice dynamics or rheology.
#
# ## Overview
#
# A slab sea ice model simplifies the full sea ice dynamics by:
# - Representing sea ice as thickness and concentration per horizontal grid point
# - Evolving thickness based on net heat fluxes and latent heat of fusion
# - Having no sea ice dynamics (no velocities, no rheology)
#
# To integrate a slab sea ice into `OceanSeaIceModel`, we need to extend several methods from ClimaOcean's
# `OceanSeaIceModels.jl` and `InterfaceComputations.jl` modules as well as Oceananigans' TimeSteppers.jl module:
#
# 1. `ComponentExchanger`: Defines how to extract state from the component for flux computations
# 2. `net_fluxes`: Returns the flux fields that will be updated by the coupling system
# 3. `sea_ice_thickness`: Returns the ice thickness field
# 4. `sea_ice_concentration`: Returns the ice concentration field
# 5. `sea_ice_top_temperature`: Returns the surface temperature field
# 6. `time_step!`: Advances the component forward in time
#
# ## Define the Slab Sea Ice Type
#
# First, we define a struct to hold the slab sea ice state and its constructor.

using Oceananigans
using Oceananigans.Fields: ConstantField

struct SlabSeaIce{G, C, H, A, T, FT, FB}
    grid :: G
    clock :: C
    thickness :: H
    concentration :: A
    top_temperature :: T
    top_heat_flux :: FT
    bottom_heat_flux :: FB
end

function SlabSeaIce(grid; clock=Clock(time=0))
    thickness = CenterField(grid)
    concentration = CenterField(grid)
    top_temperature = CenterField(grid)
    top_heat_flux = CenterField(grid)
    bottom_heat_flux = CenterField(grid)
    return SlabSeaIce(grid, clock, thickness, concentration,
                      top_temperature, top_heat_flux, bottom_heat_flux)
end

Base.summary(ice::SlabSeaIce) = "SlabSeaIce"
Base.show(io::IO, ice::SlabSeaIce) = print(io, summary(ice))

# ## Extend the InterfaceComputations.jl module
#
# The `ComponentExchanger` type contains the state variables needed for flux computations.
# For sea ice, we need to expose thickness, concentration, and velocities (zero for slab model).

using ClimaOcean.OceanSeaIceModels.InterfaceComputations

function InterfaceComputations.ComponentExchanger(ice::SlabSeaIce, exchange_grid)
    h = ice.thickness
    ℵ = ice.concentration
    hc = ConstantField(0.05)  # Consolidation thickness
    u = ZeroField()    # No dynamics
    v = ZeroField()
    return InterfaceComputations.ComponentExchanger((; u, v, hc, h, ℵ), nothing)
end

# The `net_fluxes` function returns the flux fields that will be updated by the coupling system.

function InterfaceComputations.net_fluxes(ice::SlabSeaIce)
    return (top    = (; heat = ice.top_heat_flux),
            bottom = (; heat = ice.bottom_heat_flux))
end

# ## Extend the OceanSeaIceModels.jl module
#
# In the OceanSeaIceModels.jl module, we define the sea ice state accessors and properties.

using ClimaOcean.OceanSeaIceModels
using ClimaOcean.SeaIces

OceanSeaIceModels.sea_ice_thickness(ice::SlabSeaIce) = ice.thickness
OceanSeaIceModels.sea_ice_concentration(ice::SlabSeaIce) = ice.concentration
OceanSeaIceModels.sea_ice_top_temperature(ice::SlabSeaIce) = ice.top_temperature

# The `update_net_fluxes!` function computes net fluxes and applies them to previously defined `net_fluxes` containers.

OceanSeaIceModels.update_net_fluxes!(coupled_model, ice::SlabSeaIce) =
    SeaIces.update_net_sea_ice_fluxes!(coupled_model, ice, ice.grid)

# ## Extend the TimeSteppers.jl module
#
# The `time_step!` method is called by the coupled model to advance the component forward in time.
# For slab sea ice, this method evolves thickness based on heat fluxes and latent heat of fusion.
#
# The physics is:
# ```math
# \frac{\partial h}{\partial t} = -\frac{Q_{top} + Q_{bottom}}{\rho_i \cdot L_f}
# ```
# where $\rho_i$ is ice density (~917 kg/m³) and $L_f$ is latent heat of fusion (~334,000 J/kg).

import Oceananigans.TimeSteppers: time_step!
using Oceananigans.TimeSteppers: tick!

const ρᵢ = 917.0      # Ice density (kg/m³)
const Lf = 334000.0   # Latent heat of fusion (J/kg)

function time_step!(ice::SlabSeaIce, Δt)
    tick!(ice.clock, Δt)

    h = parent(ice.thickness)
    ℵ = parent(ice.concentration)
    Q_top = parent(ice.top_heat_flux)
    Q_bot = parent(ice.bottom_heat_flux)

    # Thickness change: dh/dt = -(Q_top + Q_bot) / (ρᵢ * Lf)
    # Positive Q means heat into ice → melting → negative dh
    @. h -= (Q_top + Q_bot) * Δt / (ρᵢ * Lf)

    # Ensure non-negative thickness
    @. h = max(h, 0.0)

    # If thickness goes to zero, concentration goes to zero
    @. ℵ = ifelse(h > 0.01, ℵ, 0.0)

    return nothing
end

# ## Complete Example: Coupling Slab Sea Ice with Slab Ocean and JRA55
#
# Here's a complete example showing how to use the slab sea ice in a coupled simulation.
# We use the JRA55 reanalysis for the atmosphere and the ECCO4Monthly dataset to initialize.

using ClimaOcean
using Oceananigans.Units
using Dates

arch = GPU()

# Create a 2D grid for slab models
grid = TripolarGrid(arch, size=(180, 90, 1), z=(-50, 0))
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=5)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

# We need a SlabOcean - import from slab_ocean.jl example or define here
include("slab_ocean.jl")

slab_ocean = SlabOcean(grid)
set!(slab_ocean.temperature, Metadatum(:temperature, dataset=ECCO4Monthly()))

# Create slab sea ice and initialize with ECCO data
slab_ice = SlabSeaIce(grid)
set!(slab_ice.thickness, Metadatum(:sea_ice_thickness, dataset=ECCO4Monthly()))
set!(slab_ice.concentration, Metadatum(:sea_ice_concentration, dataset=ECCO4Monthly()))

# Load atmosphere
atmosphere = JRA55PrescribedAtmosphere(arch)

# Create coupled model - the framework handles all flux computations!
coupled_model = OceanSeaIceModel(slab_ocean, slab_ice; atmosphere)

# Run simulation
simulation = Simulation(coupled_model, Δt=1hour, stop_time=30days)
run!(simulation)

# ## Visualization

using CairoMakie

T = interior(slab_ocean.temperature, :, :, 1)
h = interior(slab_ice.thickness, :, :, 1)
ℵ = interior(slab_ice.concentration, :, :, 1)

fig = Figure(size=(1200, 400), fontsize=16)

ax1 = Axis(fig[1, 1], title="Slab Ocean Temperature (°C)")
ax2 = Axis(fig[1, 2], title="Slab Sea Ice Thickness (m)")
ax3 = Axis(fig[1, 3], title="Slab Sea Ice Concentration")

hm1 = heatmap!(ax1, Array(T); colormap=:thermal, colorrange=(-2, 30), nan_color=:lightgray)
hm2 = heatmap!(ax2, Array(h); colormap=:ice, colorrange=(0, 3), nan_color=:lightgray)
hm3 = heatmap!(ax3, Array(ℵ); colormap=:deep, colorrange=(0, 1), nan_color=:lightgray)

Colorbar(fig[2, 1], hm1, vertical=false)
Colorbar(fig[2, 2], hm2, vertical=false)
Colorbar(fig[2, 3], hm3, vertical=false)

for ax in [ax1, ax2, ax3]
    hidedecorations!(ax)
end

save("slab_sea_ice.png", fig)
nothing #hide

# ![](slab_sea_ice.png)

# ## Summary
#
# The slab sea ice example demonstrates how a simplified component can be integrated into the coupling framework.
# The key insight is that the coupling system handles flux computations generically;
# your component just needs to provide the right interface methods:
#
# | Method | Purpose |
# |--------|---------|
# | `ComponentExchanger(component, grid)` | Define what state variables the component exposes |
# | `net_fluxes(component)` | Define where computed fluxes should be stored |
# | `time_step!(component, Δt)` | Advance the component in time |
# | `sea_ice_thickness(ice)` | Return the thickness field |
# | `sea_ice_concentration(ice)` | Return the concentration field |
# | `sea_ice_top_temperature(ice)` | Return the surface temperature field |
