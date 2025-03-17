using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf
using ClimaOcean.ECCO

arch = CPU()

Nx = 144
Ny = 60
Nz = 40

depth = 6000meters
z_faces = exponential_z_faces(; Nz, depth)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces,
                             latitude  = (-75, 75),
                             longitude = (0, 360))

ocean = ocean_simulation(grid)

date = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T=Metadata(:temperature; dates=date, dataset=ECCO4Monthly()),
                  S=Metadata(:salinity; dates=date, dataset=ECCO4Monthly()))

radiation = Radiation(arch)

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))

coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

simulation = Simulation(coupled_model; Δt=30minutes, stop_time=10days)

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

    msg = @sprintf("Iter: %d, simulation time: %s, atmosphere time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(atmosphere.clock.time), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    @info msg

    wall_time[] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, TimeInterval(1days))

outputs = merge(ocean.model.tracers, ocean.model.velocities)

avg = []
avg_z = []
int = []

ρₒ = simulation.model.interfaces.ocean_properties.reference_density
cₚ = simulation.model.interfaces.ocean_properties.heat_capacity
S₀ = 35 #g/kg

for key in keys(ocean.model.tracers)
    push!(avg, (Average(ocean.model.tracers[key])))
    push!(avg_z, (Average(ocean.model.tracers[key]), dims = (1,2)))
    push!(int, (Integral(ocean.model.tracers[key]), dims = (1,2)))
end

for key in keys(ocean.model.velocities)
    push!(avg_z, (Average(ocean.model.velocities[key]), dims = (1,2)))
end

c = CenterField(grid)
volmask =  set!(c, 1)

dv = (Integral(volmask))

# Globally Averaged tracers
avg_tracer_outputs = merge((; T_avg=avg[1], S_avg=avg[2], e_avg=avg[3]))
# Depth-integrated tracers
int_tracer_outputs = merge((; OHC=ρₒ*cₚ*int[1], S_int=(ρₒ/S₀)*int[2], e_int=int[3], vol_int = dv))

depth_avg_tracer_outputs = merge((; T_avgz=avg_z[1], S_avgz=avg_z[2], e_avgz=avg_z[3]))
depth_avg_velocity_outputs = merge((; u_avgz=avg_z[4], v_avgz=avg_z[5], w_avgz=avg_z[6]))
depth_avg_outputs = merge(depth_avg_tracer_outputs, depth_avg_velocity_outputs)


simulation.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs;
                                                  schedule = TimeInterval(5days),
                                                  filename = "global_surface_fields",
                                                  indices = (:, :, grid.Nz),
                                                  with_halos = true,
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32})

simulation.output_writers[:global_avg] = JLD2OutputWriter(ocean.model, avg_tracer_outputs;
                                                  schedule = TimeInterval(1days),
                                                  filename = "averaged_tracer_data",
                                                  overwrite_existing = true)

simulation.output_writers[:global_depth_avg] = JLD2OutputWriter(ocean.model, depth_avg_tracer_outputs;
                                                  schedule = TimeInterval(1days),
                                                  filename = "depth_averaged_tracer_data",
                                                  overwrite_existing = true)


simulation.output_writers[:global_int] = JLD2OutputWriter(ocean.model, int_tracer_outputs;
                                                  schedule = TimeInterval(1days),
                                                  filename = "integrated_tripolar_data",
                                                  overwrite_existing = true)

run!(simulation)