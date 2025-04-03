using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using ClimaOcean.OceanSimulations
using ClimaOcean.ECCO
using ClimaOcean.DataWrangling
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf

using CUDA
CUDA.device!(1)
arch = GPU()

r_faces = ClimaOcean.exponential_z_faces(; Nz=100, h=30, depth=6200)
z_faces = MutableVerticalDiscretization(r_faces)

Nx = 1440 # longitudinal direction 
Ny = 700  # meridional direction 
Nz = length(r_faces) - 1

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

#####
##### A Propgnostic Ocean model
#####

# A very diffusive ocean
momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; substeps=70)
closure = ClimaOcean.OceanSimulations.default_ocean_closure()

ocean = ocean_simulation(grid;
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

# Remember to pass the SSS as a bottom bc to the sea ice!
SSS = view(ocean.model.tracers.S.data, :, :, grid.Nz)
bottom_heat_boundary_condition = IceWaterThermalEquilibrium(SSS)

# Default dynamics
dynamics = ClimaOcean.SeaIceSimulations.default_sea_ice_dynamics(grid; ocean)

sea_ice = sea_ice_simulation(grid; bottom_heat_boundary_condition, dynamics, advection=WENO(order=7))

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset),
                    ℵ=Metadatum(:sea_ice_concentration; dataset))

#####
##### A Prescribed Atmosphere model
#####

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(40))
radiation  = Radiation()

#####
##### Arctic coupled model
#####

omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
omip = Simulation(arctic, Δt=1minutes, stop_time=30days)

# Figure out the outputs....

wall_time = Ref(time_ns())

using Statistics

function progress(sim)
    sea_ice = sim.model.sea_ice
    hmax = maximum(sea_ice.model.ice_thickness)
    ℵmax = maximum(sea_ice.model.ice_concentration)
    hmean = mean(sea_ice.model.ice_thickness)
    ℵmean = mean(sea_ice.model.ice_concentration)
    Tmax = maximum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)
    Tmin = minimum(sim.model.interfaces.atmosphere_sea_ice_interface.temperature)

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg1 = @sprintf("time: %s, iteration: %d, Δt: %s, ", prettytime(sim), iteration(sim), prettytime(sim.Δt))
    msg2 = @sprintf("max(h): %.2e m, max(ℵ): %.2e ", hmax, ℵmax)
    msg3 = @sprintf("mean(h): %.2e m, mean(ℵ): %.2e ", hmean, ℵmean)
    msg4 = @sprintf("extrema(T): (%.2f, %.2f) ᵒC, ", Tmax, Tmin)
    msg5 = @sprintf("wall time: %s \n", prettytime(step_time))

    @info msg1 * msg2 * msg3 * msg4 * msg5

     wall_time[] = time_ns()

     return nothing
end

# And add it as a callback to the simulation.
add_callback!(omip, progress, IterationInterval(10))

run!(omip)

# The full OMIP cycle!
omip.stop_time = 60*365days
omip.Δt = 600

run!(omip)