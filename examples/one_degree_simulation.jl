using ClimaOcean
using ClimaOcean.ECCO: ECCO4Monthly, NearestNeighborInpainting
using OrthogonalSphericalShellGrids
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf

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
                               halo = (5, 5, 4),
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

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(tampered_bottom_height); active_cells_map=true)

#####
##### Closures
#####

# We don't need any diffusion: κ_symmetric = 0!
eddy_closure = Oceananigans.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(κ_skew=1000) 
vertical_mixing = ClimaOcean.OceanSimulations.default_ocean_closure()

closure = (eddy_closure, vertical_mixing)

#####
##### Restoring
#####

restoring_rate  = 1 / 10days
z_below_surface = z_faces[end-1]

mask = LinearlyTaperedPolarMask(southern=(-80, -70), northern=(70, 90), z=(z_below_surface, 0))

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(1993, 11, 1)
temperature = ECCOMetadata(:temperature; dates, version=ECCO4Monthly(), dir="./")
salinity    = ECCOMetadata(:salinity;    dates, version=ECCO4Monthly(), dir="./")

FT = ECCORestoring(temperature, grid; mask, rate=restoring_rate)
FS = ECCORestoring(salinity,    grid; mask, rate=restoring_rate)
forcing = (T=FT, S=FS)

#####
##### Ocean simulation
##### 

momentum_advection = WENOVectorInvariant(vorticity_order=3)
tracer_advection   = Centered()

using Oceananigans.TurbulenceClosures: IsopycnalSkewSymmetricDiffusivity,
                                       ExplicitTimeDiscretization,
                                       DiffusiveFormulation

using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity

numerical_closure = HorizontalScalarDiffusivity(ν=5e3)
eddy_closure = IsopycnalSkewSymmetricDiffusivity(κ_skew=1e3, κ_symmetric=1e3, skew_flux_formulation=DiffusiveFormulation())
vertical_mixing = CATKEVerticalDiffusivity()

closure = (eddy_closure, numerical_closure, vertical_mixing)

# Spacings still don't work correctly?
free_surface = SplitExplicitFreeSurface(grid; substeps=30)

# Should we add a side drag since this is at a coarser resolution?
ocean = ocean_simulation(grid;
                         momentum_advection,
                         tracer_advection,
                         closure,
                         forcing,
                         free_surface)

#####
##### Atmospheric forcing
#####

radiation  = Radiation(arch)
atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(20))

#####
##### Coupled simulation
#####

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation) 
simulation = Simulation(coupled_model; Δt=5minutes, stop_time=10days)

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

outputs = merge(ocean.model.tracers, ocean.model.velocities)
ocean.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                  schedule = TimeInterval(5days),
                                                  filename = "global_surface_fields",
                                                  indices = (:, :, grid.Nz),
                                                  with_halos = true,
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32})

run!(simulation)

simulation.Δt = 30minutes
simulation.stop_time = 720days

run!(simulation)

# ## A pretty movie
#
# It's time to make a pretty movie of the simulation. First we plot a snapshot:

u = FieldTimeSeries("global_surface_fields.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("global_surface_fields.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("global_surface_fields.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("global_surface_fields.jld2", "e"; backend = OnDisk())

times = u.times
Nt = length(times)

n = Observable(Nt)

land = interior(T.grid.immersed_boundary.bottom_height) .>= 0

Tn = @lift begin
    Tn = interior(T[$n])
    Tn[land] .= NaN
    view(Tn, :, :, 1)
end

en = @lift begin
    en = interior(e[$n])
    en[land] .= NaN
    view(en, :, :, 1)
end

un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    sn = interior(s)
    sn[land] .= NaN
    view(sn, :, :, 1)
end

fig = Figure(size = (800, 1200))

axs = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axT = Axis(fig[2, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axe = Axis(fig[3, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

hm = heatmap!(axs, sn, colorrange = (0, 0.5), colormap = :deep, nan_color=:lightgray)
Colorbar(fig[1, 2], hm, label = "Surface speed (m s⁻¹)")

hm = heatmap!(axT, Tn, colorrange = (-1, 30), colormap = :magma, nan_color=:lightgray)
Colorbar(fig[2, 2], hm, label = "Surface Temperature (ᵒC)")

hm = heatmap!(axe, en, colorrange = (0, 1e-3), colormap = :solar, nan_color=:lightgray)
Colorbar(fig[3, 2], hm, label = "Turbulent Kinetic Energy (m² s⁻²)")
save("snapshot.png", fig)
nothing #hide

# ![](snapshot.png)

# And now a movie:

record(fig, "near_global_ocean_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](near_global_ocean_surface.mp4)