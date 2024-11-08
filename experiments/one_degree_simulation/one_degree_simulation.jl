using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf

arch = GPU()
z = exponential_z_faces(Nz=40, depth=6000)
Nx = 360
Ny = 180
Nz = length(z_faces) - 1

grid = TripolarGrid(arch; z, size = (Nx, Ny, Nz), north_poles_latitude=55, first_pole_longitude=70)

bottom_height = regrid_bathymetry(grid;
                                  minimum_depth = 10,
                                  interpolation_passes = 5,
                                  major_basins = 3)

gm = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1000, κ_symmetric=1000)
catke = Oceananigans.TurbulenceClosures.CATKEVerticalDiffusivity()
closure = (gm, catke)

momentum_advection = VectorInvariant()
tracer_advection = Centered(order=2)
ocean = ocean_simulation(grid; momentum_advection, tracer_advection, closure)
start_date = DateTimeProlepticGregorian(1993, 1, 1)

set!(ocean.model,
     T = ECCOMetadata(:temperature; date=start_date),
     S = ECCOMetadata(:salinity; date=start_date))

atmosphere = JRA55_prescribed_atmosphere(arch; backend=JRA55NetCDFBackend(41))
radiation = Radiation(arch)
sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
simulation = Simulation(coupled_model; Δt=10minutes, stop_time=30days)

wall_time = Ref(time_ns)

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    Tmax = maximum(T)
    Tmin = minimum(T)
    umax = maximum(abs, u), maximum(abs, v), maximum(abs, w)
    step_time = 1e-9 * (time_ns() - wall_time[])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T): %.2f, min(T): %.2f, wtime: %s \n",
                   prettytime(ocean.model.clock.time),
                   ocean.model.clock.iteration,
                   prettytime(ocean.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[] = time_ns()

     return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

run!(simulation)

