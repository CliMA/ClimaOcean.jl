using ClimaOcean
using ClimaOcean.ECCO: ECCO4Monthly, NearestNeighborInpainting
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf
using CUDA: @allowscalar, device!

using Oceananigans.Grids: znode

arch = GPU()

#####
##### Grid and Bathymetry
#####

Nx = 360
Ny = 180
Nz = 100

z_faces = exponential_z_faces(; Nz, depth=5000, h=34)

underlying_grid = TripolarGrid(arch;
                               size = (Nx, Ny, Nz),
                               z = z_faces,
                               first_pole_longitude = 70,
                               north_poles_latitude = 55)

bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = 75,
                                  major_basins = 2)

# Open Gibraltar strait 
# TODO: find a better way to do this
tampered_bottom_height = deepcopy(bottom_height)
view(tampered_bottom_height, 102:103, 124, 1) .= -400

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(tampered_bottom_height))

#####
##### Closures
#####

gm = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1000, κ_symmetric=1000)
catke = ClimaOcean.OceanSimulations.default_ocean_closure()
viscous_closure = Oceananigans.TurbulenceClosures.HorizontalScalarDiffusivity(ν=2000)

closure = (gm, catke, viscous_closure)

#####
##### Restoring
#####

restoring_rate  = 1 / 2days
z_below_surface = @allowscalar znode(1, 1, grid.Nz, grid, Center(), Center(), Face())

mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(z_below_surface, 0))

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(1993, 11, 1)
temperature = ECCOMetadata(:temperature; dates, version=ECCO4Monthly(), dir="./")
salinity    = ECCOMetadata(:salinity;    dates, version=ECCO4Monthly(), dir="./")

# inpainting = NearestNeighborInpainting(30) should be enough to fill the gaps near bathymetry
FT = ECCORestoring(temperature, grid; mask, rate=restoring_rate, inpainting=NearestNeighborInpainting(50))
FS = ECCORestoring(salinity, grid;    mask, rate=restoring_rate, inpainting=NearestNeighborInpainting(50))
forcing = (T=FT, S=FS)

#####
##### Ocean simulation
##### 

momentum_advection = VectorInvariant()
tracer_advection   = Centered(order=2)

# Should we add a side drag since this is at a coarser resolution?
ocean = ocean_simulation(grid; momentum_advection, tracer_advection,
                         closure, forcing,
                         tracers = (:T, :S, :e))

set!(ocean.model, T=ECCOMetadata(:temperature; dates=first(dates)),
                  S=ECCOMetadata(:salinity;    dates=first(dates)))

#####
##### Atmospheric forcing
#####

radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(20))

#####
##### Coupled simulation
#####

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation) 
simulation = Simulation(coupled_model; Δt=15minutes, stop_time=2*365days)

#####
##### Run it!
##### 

wall_time = Ref(time_ns())

function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T
    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    @info @sprintf("Time: %s, n: %d, Δt: %s, max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s \n",
                   prettytime(sim), iteration(sim), prettytime(sim.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[] = time_ns()

     return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)
