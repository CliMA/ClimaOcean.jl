using ClimaOcean
using ClimaOcean.ECCO: ECCO4Monthly
using Oceananigans.OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Units
using Dates
using Printf

filename = "tidal_potential_jra55.jld2"
url = "https://www.dropbox.com/scl/fi/5qu30pcuy1z5q0s5xpwaj/" *
      "tidal_potential_jra55.jld2?rlkey=2597z4vkiuq3stjsry45rpfyd&st=n42stc80&dl=0"

if !isfile(filename)
    Base.download(url, filename)
end

Oceananigans.defaults.FloatType = Float64
arch = GPU()
Nx = 360 * 4
Ny = 170 * 4
Nz = 60
z_faces = exponential_z_faces(; Nz, depth=6000, h=34)
underlying_grid = TripolarGrid(arch; size=(Nx, Ny, Nz), halo=(7, 7, 7), z=z_faces)

bottom_height = regrid_bathymetry(underlying_grid; minimum_depth=20, interpolation_passes=40, major_basins=3)
#view(bottom_height, 73:78, 88:89, 1) .= -1000 # open Gibraltar strait 

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map=true)

dates = DateTime(1993, 1, 1) : Month(1) : DateTime(1993, 11, 1)
mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(-100, 0))
temperature = ECCOMetadata(:temperature; dates, version=ECCO4Monthly())
salinity    = ECCOMetadata(:salinity;    dates, version=ECCO4Monthly())
rate = 1/10days
FT = ECCORestoring(temperature, grid; mask, rate)
FS = ECCORestoring(salinity, grid; mask, rate)
forcing = (T=FT, S=FS)

#=
catke = ClimaOcean.OceanSimulations.default_ocean_closure()
viscous_closure = Oceananigans.TurbulenceClosures.HorizontalScalarDiffusivity(ν=2000)
gm = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1000, κ_symmetric=1000)
closure = (gm, catke, viscous_closure)
momentum_advection = VectorInvariant()
tracer_advection = Centered(order=2)
free_surface = SplitExplicitFreeSurface(grid; substeps=70)
ocean = ocean_simulation(grid; momentum_advection, tracer_advection, free_surface, forcing)
=#

free_surface = SplitExplicitFreeSurface(grid; substeps=70)
ocean = ocean_simulation(grid; forcing, free_surface, bottom_drag_coefficient=0.01)
                         
set!(ocean.model, T=ECCOMetadata(:temperature; dates=first(dates)),
                  S=ECCOMetadata(:salinity;    dates=first(dates)))

radiation  = Radiation(arch)
tidal_potential = FieldTimeSeries("tidal_potential_jra55.jld2", "Φ"; architecture=GPU(), backend=InMemory(41))
boundary_conditions = FieldBoundaryConditions(tidal_potential.grid, (Center, Center, Nothing))

tidal_potential = FieldTimeSeries("tidal_potential_jra55.jld2", "Φ"; architecture=GPU(), backend=InMemory(41), boundary_conditions) 
Oceananigans.BoundaryConditions.fill_halo_regions!(tidal_potential)
atmosphere = JRA55PrescribedAtmosphere(arch; tidal_potential, backend=JRA55NetCDFBackend(41))

# atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))
coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation) 
simulation = Simulation(coupled_model; Δt=1minutes, stop_time=5days)

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

u, v, w = ocean.model.velocities
ζ = ∂x(v) - ∂y(u)
s = @at (Center, Center, Center) sqrt(u^2 + v^2)
outputs = merge(ocean.model.tracers, ocean.model.velocities, (; ζ, s))
simulation.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                       schedule = TimeInterval(3hours),
                                                       filename = "tidally_forced_surface_fields",
                                                       indices = (:, :, grid.Nz),
                                                       overwrite_existing = true,
                                                       array_type = Array{Float32})

simulation.output_writers[:middepth] = JLD2OutputWriter(ocean.model, outputs;
                                                        schedule = TimeInterval(3hours),
                                                        filename = "tidally_forced_deep_fields",
                                                        indices = (:, :, Int(grid.Nz/2)),
                                                        overwrite_existing = true,
                                                        array_type = Array{Float32})

run!(simulation)

