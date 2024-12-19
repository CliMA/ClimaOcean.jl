using Printf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: φnode
using ClimaOcean
using CairoMakie

using CFTime
using Dates

using ClimaOcean.ECCO

#things to do
# 1. git pull
# 2. update library
# 3. add CairoMakie
# 4. add oceananigans#main
# 5. ECCO_USERNAME=francispoulin ECCO_PASSWORD=???????????? julia --project

# 1) No output, need to add
# 2) Restoring force (work in progress)

arch = GPU() 

z_faces = exponential_z_faces(Nz=40, depth=6000)

Nx = 1440
Ny = 600
Nz = length(z_faces) - 1

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces, 
                             latitude  = (-80, -20),
                             longitude = (0, 360))

bottom_height = regrid_bathymetry(grid; 
                                  minimum_depth = 10,
                                  interpolation_passes = 5,
                                  major_basins = 1)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height), active_cells_map=true) 

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(1993, 5, 1)

#
# Restoring force 
#               φS                   φN             -20
# -------------- | ------------------ | ------------ |
# no restoring   0    linear mask     1   mask = 1   1
#

const φN₁ = -22
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

Fu = Forcing(u_restoring; discrete_form=true, parameters=(; rate=1/2days))
Fv = Forcing(v_restoring; discrete_form=true, parameters=(; rate=1/2days))

FT = ECCORestoring(:temperature, grid; dates, rate=1/2days, mask=tracer_mask)
FS = ECCORestoring(:salinity,    grid; dates, rate=1/2days, mask=tracer_mask)

forcing = (T=FT, S=FS, u=Fu, v=Fv)

ocean = ocean_simulation(grid; forcing)
model = ocean.model

set!(model, 
     T = ECCOMetadata(:temperature; dates = dates[1]),
     S = ECCOMetadata(:salinity;    dates = dates[1]))
     
backend    = JRA55NetCDFBackend(41) 
atmosphere = JRA55PrescribedAtmosphere(arch; backend)
radiation  = Radiation()

coupled_model      = OceanSeaIceModel(ocean; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model; Δt=10, stop_time = 10days)

wall_time = [time_ns()]

function progress(sim) 
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T): %.2f, min(T): %.2f, wtime: %s \n",
                   prettytime(ocean.model.clock.time),
                   ocean.model.clock.iteration,
                   prettytime(ocean.Δt),
                   umax..., Tmax, Tmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))

ocean.output_writers[:surface] = JLD2OutputWriter(model, merge(model.tracers, model.velocities);
                                                  schedule = TimeInterval(1days),
                                                  filename = "surface",
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
# integration with a maximum time step of 1.5 minutes should be sufficient to dissipate spurious
# initialization shocks.
# We use an adaptive time step that maintains the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) equal to 0.1.
# For this scope, we use the Oceananigans utility `conjure_time_step_wizard!` (see Oceanigans's documentation).

ocean.stop_time = 10days
conjure_time_step_wizard!(ocean; cfl = 0.1, max_Δt = 90, max_change = 1.1)
run!(coupled_simulation)
nothing #hide

# ### Running the simulation
#
# Now that the simulation has spun up, we can run it for the full 100 days.
# We increase the maximum time step size to 10 minutes and let the simulation run for 100 days.
# This time, we set the CFL in the time_step_wizard to be 0.25 as this is the maximum recommended CFL to be
# used in conjunction with Oceananigans' hydrostatic time-stepping algorithm ([two step Adams-Bashfort](https://en.wikipedia.org/wiki/Linear_multistep_method))

ocean.stop_time = 2*365days
coupled_simulation.stop_time = 100days
conjure_time_step_wizard!(ocean; cfl = 0.25, max_Δt = 10minutes, max_change = 1.1)
run!(coupled_simulation)
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

iter = Observable(Nt)

Ti = @lift begin
     Ti = interior(T[$iter], :, :, 1)
     Ti[Ti .== 0] .= NaN
     Ti
end

ei = @lift begin
     ei = interior(e[$iter], :, :, 1)
     ei[ei .== 0] .= NaN
     ei
end

si = @lift begin
     s = Field(sqrt(u[$iter]^2 + v[$iter]^2))
     compute!(s)
     s = interior(s, :, :, 1)
     s[s .== 0] .= NaN
     s
end


fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, si, colorrange = (0, 0.5), colormap = :deep)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Surface speed (ms⁻¹)")
hidedecorations!(ax)

CairoMakie.record(fig, "near_global_ocean_surface_s.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
end
nothing #hide
 
 # ![](near_global_ocean_surface_s.mp4)
 
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Ti, colorrange = (-1, 30), colormap = :magma)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Surface Temperature (Cᵒ)")
hidedecorations!(ax)

CairoMakie.record(fig, "near_global_ocean_surface_T.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
end
nothing #hide
 
# ![](near_global_ocean_surface_T.mp4)

fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, ei, colorrange = (0, 1e-3), colormap = :solar)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Turbulent Kinetic Energy (m²s⁻²)")
hidedecorations!(ax)

CairoMakie.record(fig, "near_global_ocean_surface_e.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
end
nothing #hide

# ![](near_global_ocean_surface_e.mp4)
