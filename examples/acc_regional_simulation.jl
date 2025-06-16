using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnode
using ClimaOcean
using ClimaOcean.ECCO

using Printf
using CairoMakie
using CFTime
using Dates

arch = GPU() 

Nx = 1440
Ny = 400
Nz = 40  

z_faces = exponential_z_faces(; Nz, depth=6000)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces, 
                             latitude  = (-80, -20),
                             longitude = (0, 360))

bottom_height = regrid_bathymetry(grid; 
                                  minimum_depth = 10,
                                  interpolation_passes = 7,
                                  major_basins = 1)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height), active_cells_map=true) 

start_date = DateTime(1993, 1, 1)
end_date   = DateTime(1993, 12, 1) 

#
# Restoring force 
#               φS                   φN             -20
# -------------- | ------------------ | ------------ |
# no restoring   0    linear mask     1   mask = 1   1
#

const φN₁ = -23
const φN₂ = -25
const φS₁ = -78
const φS₂ = -75

@inline northern_mask(φ)    = min(max((φ - φN₂) / (φN₁ - φN₂), zero(φ)), one(φ))
@inline southern_mask(φ, z) = ifelse(z > -20, 
                                     min(max((φ - φS₂) / (φS₁ - φS₂), zero(φ)), one(φ)), 
                                     zero(φ))

@inline function tracer_mask(λ, φ, z, t)
     n = northern_mask(φ)
     s = southern_mask(φ, z)
     return max(s, n)
end

@inline function u_restoring(i, j, k, grid, clock, fields, p)
     φ = φnode(i, j, k, grid, Face(), Center(), Center())
     return - p.rate * fields.u[i, j, k] * northern_mask(φ)
end

@inline function v_restoring(i, j, k, grid, clock, fields, p)
     φ = φnode(i, j, k, grid, Center(), Face(), Center())
     return - p.rate * fields.v[i, j, k] * northern_mask(φ)
end

T_meta = Metadata(:temperature; start_date, end_date, dataset = ECCO4Monthly())
S_meta = Metadata(:temperature; start_date, end_date, dataset = ECCO4Monthly())

forcing = (T = DatasetRestoring(T_meta, grid; rate=1/5days, mask=tracer_mask),
	       S = DatasetRestoring(S_meta, grid; rate=1/5days, mask=tracer_mask),
	       u = Forcing(u_restoring; discrete_form=true, parameters=(; rate=1/5days)),
	       v = Forcing(v_restoring; discrete_form=true, parameters=(; rate=1/5days)))

momentum_advection = WENOVectorInvariant()
tracer_advection   = WENO(order=7)

ocean      = ocean_simulation(grid; forcing, momentum_advection, tracer_advection)
set!(ocean.model, T=0, S=0)
backend    = JRA55NetCDFBackend(41) 
atmosphere = JRA55PrescribedAtmosphere(arch; backend)
radiation  = Radiation()
model      = ocean.model 

coupled_model = OceanSeaIceModel(ocean; atmosphere=nothing, radiation)
simulation    = Simulation(coupled_model; Δt=2minutes, stop_time = 10days)

wall_time = [time_ns()]

function progress(sim) 
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax, Tmin = maximum(T), minimum(T)
    umax = maximum(abs, u), maximum(abs, v), maximum(abs, w)
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T): %.2f, min(T): %.2f, wtime: %s \n",
                   prettytime(ocean.model.clock.time),
                   ocean.model.clock.iteration,
                   prettytime(ocean.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, TimeInterval(4hours)) 

ocean.output_writers[:surface] = JLD2Writer(model, merge(model.tracers, model.velocities);
                                            schedule = TimeInterval(5days),
                                            filename = "acc_surface_fields",
                                            indices = (:, :, grid.Nz),
                                            overwrite_existing = true,
                                            array_type = Array{Float32})
nothing #hide

# ### Spinning up the simulation
#
# As an initial condition, we have interpolated ECCO tracer fields onto our custom grid.
# The bathymetry of the original ECCO data may differ from our grid, so the initialization of the velocity
# field might cause shocks if a large time step is used.
#
# Therefore, we spin up the simulation with a small time step to ensure that the interpolated initial
# conditions adapt to the model numerics and parameterization without causing instability. A 10-day
# integration with a maximum time step of 1 minute should be sufficient to dissipate spurious
# initialization shocks.

run!(simulation)
nothing #hide

# ### Running the simulation
#
# Now that the simulation has spun up, we can run it for the full 2 years.
# We increase the maximum time step size to 10 minutes and let the simulation run for 2 years.

simulation.stop_time = 2*365days
simulation.Δt = 10minutes
run!(simulation)
nothing #hide

# ## Visualizing the results
# 
# The simulation has finished, let's visualize the results.
# In this section we pull up the saved data and create visualizations using the CairoMakie.jl package.
# In particular, we generate an animation of the evolution of surface fields:
# surface speed (s), surface temperature (T), and turbulent kinetic energy (e).

u = FieldTimeSeries("surface.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("surface.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("surface.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("surface.jld2", "e"; backend = OnDisk())

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

s = @at (Center, Center, Nothing) sqrt(un^2 + vn^2) 
s = Field(s)

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    sn = interior(s)
    sn[land] .= NaN
    view(sn, :, :, 1)
end

title = @lift string("ACC regional ocean simulation after ",
                     prettytime(times[$n] - times[1]))

λ, φ, _ = nodes(T)

fig = Figure(size = (1000, 1500))

axs = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axT = Axis(fig[2, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axe = Axis(fig[3, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

hm = heatmap!(axs, λ, φ, sn, colorrange = (0, 0.5), colormap = :deep, nan_color=:lightgray)
Colorbar(fig[1, 2], hm, label = "Surface Speed (m s⁻¹)")

hm = heatmap!(axT, λ, φ, Tn, colorrange = (-1, 30), colormap = :magma, nan_color=:lightgray)
Colorbar(fig[2, 2], hm, label = "Surface Temperature (ᵒC)")

hm = heatmap!(axe, λ, φ, en, colorrange = (0, 1e-3), colormap = :solar, nan_color=:lightgray)
Colorbar(fig[3, 2], hm, label = "Turbulent Kinetic Energy (m² s⁻²)")

Label(fig[0, :], title)

save("snapshot_acc.png", fig)
nothing #hide

# ![](snapshot.png)

# And now we make a movie:

record(fig, "acc_region_surface.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](near_global_ocean_surface.mp4)
