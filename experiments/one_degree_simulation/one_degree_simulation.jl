using ClimaOcean
using ClimaOcean.ECCO: ECCO4Monthly
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Units
using Dates
using Printf
using GLMakie

Oceananigans.defaults.FloatType = Float64
arch = CPU()
Nx = 256
Ny = 128
Nz = 32
z_faces = exponential_z_faces(; Nz, depth=6000, h=34)
underlying_grid = TripolarGrid(arch; size=(Nx, Ny, Nz), z=z_faces)

bottom_height = regrid_bathymetry(underlying_grid; minimum_depth=30, interpolation_passes=20, major_basins=1)
view(bottom_height, 73:78, 88:89, 1) .= -1000 # open Gibraltar strait 

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map=true)

catke = ClimaOcean.OceanSimulations.default_ocean_closure()
viscous_closure = Oceananigans.TurbulenceClosures.HorizontalScalarDiffusivity(ν=2000)
closure = (catke, viscous_closure)

#gm = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1000, κ_symmetric=1000)
#closure = (gm, catke, viscous_closure)

dates = DateTime(1993, 1, 1) : Month(1) : DateTime(1993, 11, 1)
mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(-100, 0))
temperature = ECCOMetadata(:temperature; dates, version=ECCO4Monthly())
salinity    = ECCOMetadata(:salinity;    dates, version=ECCO4Monthly())
rate = 1/10days
FT = ECCORestoring(temperature, grid; mask, rate)
FS = ECCORestoring(salinity, grid; mask, rate)
forcing = (T=FT, S=FS)

momentum_advection = VectorInvariant()
tracer_advection = Centered(order=2)
free_surface = SplitExplicitFreeSurface(grid; substeps=70)
ocean = ocean_simulation(grid; momentum_advection, tracer_advection, free_surface, forcing)
                         
set!(ocean.model, T=ECCOMetadata(:temperature; dates=first(dates)),
                  S=ECCOMetadata(:salinity;    dates=first(dates)))

radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))
coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation) 

#=
simulation = Simulation(coupled_model; Δt=20minutes, stop_iteration=100)

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
=#

