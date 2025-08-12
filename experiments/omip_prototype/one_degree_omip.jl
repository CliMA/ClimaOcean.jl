using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf
using Dates
using CUDA

import Oceananigans.OutputWriters: checkpointer_address

function synch!(clock1::Clock, clock2)
    # Synchronize the clocks
    clock1.time = clock2.time
    clock1.iteration = clock2.iteration
    clock1.last_Δt = clock2.last_Δt
end

synch!(model1, model2) = synch!(model1.clock, model2.clock)

arch = GPU()

Nx = 360 # longitudinal direction 
Ny = 180 # meridional direction 
Nz = 60

r_faces = ClimaOcean.ExponentialCoordinate(Nz, -6000)
z_faces = MutableVerticalDiscretization(r_faces)

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=15)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

#####
##### A Propgnostic Ocean model
#####

using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation

momentum_advection = WENOVectorInvariant(order=5)
tracer_advection   = WENO(order=5)

free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt=20minutes)

mixing_length = CATKEMixingLength(Cᵇ=0.01)
turbulent_kinetic_energy_equation = CATKEEquation(Cᵂϵ=1.0)

catke_closure = CATKEVerticalDiffusivity(; mixing_length, turbulent_kinetic_energy_equation) 
closure = (catke_closure, VerticalScalarDiffusivity(κ=1e-5, ν=1e-5))

ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         free_surface,
                         closure)

dataset = ECCO4Monthly()

set!(ocean.model, T=Metadatum(:temperature; dataset),
                  S=Metadatum(:salinity;    dataset))

#####
##### A Prognostic Sea-ice model
#####

# Default sea-ice dynamics and salinity coupling are included in the defaults
sea_ice = sea_ice_simulation(grid; advection=WENO(order=7))

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset),
                    ℵ=Metadatum(:sea_ice_concentration; dataset))

#####
##### A Prescribed Atmosphere model
#####

dir = "./forcing_data"
dataset = MultiYearJRA55()
backend = JRA55NetCDFBackend(100)

atmosphere = JRA55PrescribedAtmosphere(arch; dir, dataset, backend, include_rivers_and_icebergs=true)
radiation  = Radiation()

#####
##### An ocean-sea ice coupled model
#####
 
omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
omip = Simulation(omip, Δt=30, stop_time=60days)

# Figure out the outputs....

checkpointer_address(::SeaIceModel) = "SeaIceModel"

ocean.output_writers[:checkpointer] = Checkpointer(ocean.model,
                                                  schedule = IterationInterval(10000),
                                                  prefix = "ocean_checkpoint_onedegree",
                                                  overwrite_existing = true)

sea_ice.output_writers[:checkpointer] = Checkpointer(sea_ice.model,
                                                     schedule = IterationInterval(10000),
                                                     prefix = "sea_ice_checkpoint_onedegree",
                                                     overwrite_existing = true)

wall_time = Ref(time_ns())

using Statistics

function progress(sim)
    sea_ice = sim.model.sea_ice
    ocean   = sim.model.ocean
    hmax = maximum(sea_ice.model.ice_thickness)
    ℵmax = maximum(sea_ice.model.ice_concentration)
    Tmax = maximum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    Tmin = minimum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    umax = maximum(ocean.model.velocities.u)
    vmax = maximum(ocean.model.velocities.v)
    wmax = maximum(ocean.model.velocities.w)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    msg4 = @sprintf("extrema(T): (%.2f, %.2f) ᵒC, ", Tmax, Tmin)
    msg5 = @sprintf("maximum(u): (%.2f, %.2f, %.2f) m/s, ", umax, vmax, wmax)
    msg6 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg4 * msg5 * msg6

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(omip, progress, IterationInterval(50))

run!(omip)

omip.Δt = 15minutes
omip.stop_time = 58 * 365days

run!(omip)
