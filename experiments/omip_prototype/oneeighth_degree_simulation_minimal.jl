using ClimaOcean
using ClimaSeaIce
using Oceananigans
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids
using Oceananigans.Architectures: on_architecture
using Oceananigans.DistributedComputations: Equal
using Oceananigans.BoundaryConditions: fill_halo_regions!
using ClimaOcean.OceanSimulations
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: NearestNeighborInpainting
using ClimaSeaIce.SeaIceThermodynamics: IceWaterThermalEquilibrium
using Printf
using Dates
using CUDA
using PythonCall
using Oceananigans.BuoyancyFormulations: buoyancy, buoyancy_frequency

import Oceananigans.OutputWriters: checkpointer_address

# arch = GPU()
# arch = Distributed(GPU(), partition=Partition(1, 4), synchronized_communication=true)
arch = Distributed(GPU(); partition = Partition(y = Equal()), synchronized_communication=true)
# arch = Distributed(CPU(), partition=Partition(1, 2), synchronized_communication=true)

@info "Architecture $(arch)"

Nx = 2880 # longitudinal direction 
Ny = 1440 # meridional direction 
Nz = 100

# z_faces = ExponentialCoordinate(Nz, -6000, 0)
z_faces = ExponentialDiscretization(Nz, -6000, 0)

const z_surf = z_faces(Nz)

@info "Building grid..."
grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (7, 7, 7))

@info "Regridding bathymetry..."
bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=10)

fitted_bottom = GridFittedBottom(bottom_height)

@info "Building immersed boundary grid..."
grid = ImmersedBoundaryGrid(grid, fitted_bottom; active_cells_map=true)

@info grid
@info "Created ImmersedBoundaryGrid"

#####
##### A Propgnostic Ocean model
#####

using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation
using Oceananigans.TurbulenceClosures: RiBasedVerticalDiffusivity

momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

# free_surface = SplitExplicitFreeSurface(grid; cfl=0.8, fixed_Δt=12minutes)
free_surface = SplitExplicitFreeSurface(grid; substeps = 70)
@info "Free surface", free_surface

obl_closure = ClimaOcean.OceanSimulations.default_ocean_closure() # CATKE
closure = (obl_closure, VerticalScalarDiffusivity(κ=1e-5, ν=1e-4))

glorys_dir = joinpath(homedir(), "GLORYS_data")
mkpath(glorys_dir)

glorys_dataset = GLORYSMonthly()

@info "Building ocean component..."
ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         timestepper = :SplitRungeKutta3,
                         free_surface,
                         closure)

start_date = DateTime(1993, 1, 1)
end_date   = DateTime(2003, 4, 1)
simulation_period = Dates.value(Second(end_date - start_date))

inpainting = NearestNeighborInpainting(50)
@info "Setting initial conditions..."

Tᵢ = Metadatum(:temperature; dataset=glorys_dataset, date=start_date, dir=glorys_dir)
Sᵢ = Metadatum(:salinity;    dataset=glorys_dataset, date=start_date, dir=glorys_dir)

set!(ocean.model.tracers.T, Tᵢ; inpainting)
set!(ocean.model.tracers.S, Sᵢ; inpainting)

fill_halo_regions!(ocean.model.tracers.T)
fill_halo_regions!(ocean.model.tracers.S)

# set!(ocean.model, T=Metadatum(:temperature; dataset=glorys_dataset, date=start_date, dir=glorys_dir),
                #   S=Metadatum(:salinity;    dataset=glorys_dataset, date=start_date, dir=glorys_dir); inpainting)

@info ocean.model.clock

#####
##### A Prognostic Sea-ice model
#####

@info "Building sea-ice component..."
sea_ice = sea_ice_simulation(grid, ocean; dynamics=nothing)

@info "Setting sea-ice initial conditions..."
set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset=glorys_dataset, dir=glorys_dir),
                    ℵ=Metadatum(:sea_ice_concentration; dataset=glorys_dataset, dir=glorys_dir), inpainting = nothing)

#####
##### A Prescribed Atmosphere model
#####

jra55_dir = joinpath(homedir(), "JRA55_data")
mkpath(jra55_dir)
dataset = MultiYearJRA55()
jra55_backend = JRA55NetCDFBackend(10)

@info "Building atmospheric forcing..."
atmosphere = JRA55PrescribedAtmosphere(arch; dir=jra55_dir, dataset=jra55_backend, backend=jra55_backend, include_rivers_and_icebergs=true, start_date)
radiation  = Radiation()

#####
##### An ocean-sea ice coupled model
#####

@info "Building coupled ocean-sea ice model..."
omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
omip = Simulation(omip, Δt=10minutes, stop_time=60days)

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

add_callback!(omip, progress, IterationInterval(1))

@info "Starting simulation..."
run!(omip)

omip.Δt = 10minutes
omip.stop_time = simulation_period

run!(omip)