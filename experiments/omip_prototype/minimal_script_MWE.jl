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
using Oceananigans.TurbulenceClosures: DiffusiveFormulation, AdvectiveFormulation
using Oceananigans.BuoyancyFormulations: buoyancy, buoyancy_frequency

import Oceananigans.OutputWriters: checkpointer_address

# arch = GPU()
# arch = Distributed(GPU(), partition=Partition(1, 4), synchronized_communication=true)
arch = Distributed(CPU(), partition=Partition(1, 4), synchronized_communication=true)
@info "Architecture $(arch)"

Nx = 2880 # longitudinal direction 
Ny = 1440 # meridional direction 
Nz = 100

z_faces = ExponentialDiscretization(Nz, -6000, 0)

const z_surf = z_faces(Nz)
halo_size = 7

grid = TripolarGrid(arch;
                    size = (Nx, Ny, Nz),
                    z = z_faces,
                    halo = (halo_size, halo_size, halo_size))

@info "halo size: $(halo_size)"

bottom_height = regrid_bathymetry(grid; minimum_depth=15, major_basins=1, interpolation_passes=10)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; substeps = 70)

obl_closure = ClimaOcean.OceanSimulations.default_ocean_closure()
closure = (obl_closure, VerticalScalarDiffusivity(κ=1e-5, ν=1e-4))

dir = joinpath(homedir(), "forcing_data_1deg_minimal_multi40_backend3")
mkpath(dir)

dataset = EN4Monthly()
start_date = DateTime(1993, 1, 1)
end_date = start_date + Month(3)

@inline mask(x, y, z, t) = z ≥ z_surf - 1
Smetadata = Metadata(:salinity; dataset, dir)

ocean = ocean_simulation(grid; Δt=1minutes,
                         momentum_advection,
                         tracer_advection,
                         timestepper = :SplitRungeKutta3,
                         free_surface,
                         closure)

dataset = EN4Monthly()

set!(ocean.model, T=Metadatum(:temperature; dataset, start_date, dir),
                  S=Metadatum(:salinity;    dataset, start_date, dir))

@info ocean.model.clock

sea_ice = sea_ice_simulation(grid, ocean; dynamics=nothing)

set!(sea_ice.model, h=Metadatum(:sea_ice_thickness;     dataset=ECCO4Monthly(), dir, date=start_date),
                    ℵ=Metadatum(:sea_ice_concentration; dataset=ECCO4Monthly(), dir, date=start_date))

jra55_dir = joinpath(homedir(), "JRA55_data")
dataset = MultiYearJRA55()
backend = JRA55NetCDFBackend(3)

@info "dataset: $dataset"

atmosphere = JRA55PrescribedAtmosphere(arch; dir=jra55_dir, dataset, backend, include_rivers_and_icebergs=true, start_date, end_date)

@info "atmosphere: $atmosphere"
radiation  = Radiation()

@info "Setting up Ocean-SeaIce-Atmosphere model..."
omip = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

@info "Setting up OMIP simulation..."
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

# And add it as a callback to the simulation.
add_callback!(omip, progress, IterationInterval(1))

@info "Starting OMIP simulation..."
run!(omip)