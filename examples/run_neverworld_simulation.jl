using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.ImmersedBoundaries: PartialCellBottom, GridFittedBottom
using ClimaOcean.IdealizedSimulations: neverworld_simulation
using ClimaOcean.VerticalGrids: stretched_vertical_faces, PowerLawStretching
using Printf
using CUDA

closure = CATKEVerticalDiffusivity(minimum_turbulent_kinetic_energy = 1e-6,
                                   minimum_convective_buoyancy_flux = 1e-11)

z = stretched_vertical_faces(surface_layer_Δz = 8,
                             surface_layer_height = 32,
                             stretching = PowerLawStretching(1.02),
                             maximum_Δz = 400.0,
                             minimum_depth = 4000)

simulation = neverworld_simulation(GPU(); z,
                                   horizontal_resolution = 1/8,
                                   longitude = (0, 60),
                                   latitude = (-70, 0),
                                   time_step = 5minutes,
                                   stop_time = 4 * 360days,
                                   closure)

model = simulation.model
grid = model.grid

@show grid

start_time = Ref(time_ns())
previous_model_time = Ref(time(simulation))

function progress(sim) 
    b = sim.model.tracers.b
    e = sim.model.tracers.e
    u, v, w = sim.model.velocities

    msg = @sprintf("Iter: %d, time: %s, extrema(b): (%6.2e, %6.2e)",
                   iteration(sim), prettytime(sim), minimum(b), maximum(b))

    msg *= @sprintf(", max(e): %6.2e", maximum(e))

    msg *= @sprintf(", max|u|: %6.2e, max|w|: %6.2e",
                    maximum(maximum(abs, q) for q in (u, v, w)), maximum(abs, w))

    try 
        κᶜ = sim.model.diffusivity_fields.κᶜ
        msg *= @sprintf(", max(κᶜ): %6.2e", maximum(κᶜ))
    catch
    end

    elapsed = 1e-9 * (time_ns() - start_time[])
    elapsed_model_time = time(sim) - previous_model_time[]
    SYPD = (elapsed_model_time/360days) / (elapsed/day)

    msg *= @sprintf(", wall time: %s, SYPD: %.1f", prettytime(elapsed), SYPD)
    start_time[] = time_ns()
    previous_model_time[] = time(sim)

    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Set up output
Nx, Ny, Nz = size(grid)
Δt = simulation.Δt
Δt_minutes = round(Int, Δt / minutes) 
ib_str = grid.immersed_boundary isa PartialCellBottom ? "partial_cells" : "full_cells"
output_suffix = "$(Nx)_$(Ny)_$(Nz)_dt$(Δt_minutes)_$(ib_str).jld2"
fine_output_frequency = 1day
i = round(Int, Nx/10) # index for yz-sliced output

z = znodes(grid, Face(), with_halos=true)

K = CUDA.@allowscalar [Nz,
                       searchsortedfirst(z, -100),
                       searchsortedfirst(z, -400)]

i = round(Int, Nx/10) # index for yz-sliced output

diffusivity_fields = (; κᶜ = model.diffusivity_fields.κᶜ)
outputs = merge(model.velocities, model.tracers, diffusivity_fields)
zonally_averaged_outputs = NamedTuple(n => Average(outputs[n], dims=1) for n in keys(outputs))

simulation.output_writers[:yz] = JLD2OutputWriter(model, outputs;
                                                  schedule = TimeInterval(fine_output_frequency),
                                                  filename = "neverworld_yz_" * output_suffix,
                                                  indices = (i, :, :),
                                                  with_halos = true,
                                                  overwrite_existing = true)

simulation.output_writers[:zonal] = JLD2OutputWriter(model, zonally_averaged_outputs;
                                                     schedule = TimeInterval(fine_output_frequency),
                                                     filename = "neverworld_zonal_average_" * output_suffix,
                                                     with_halos = true,
                                                     overwrite_existing = true)

for (n, k) in enumerate(K)
    name = Symbol(:xy, n)
    simulation.output_writers[name] = JLD2OutputWriter(model, outputs;
                                                       schedule = TimeInterval(fine_output_frequency),
                                                       filename = "neverworld_xy$(n)_" * output_suffix,
                                                       indices = (:, :, Nz),
                                                       with_halos = true,
                                                       overwrite_existing = true)
end

simulation.output_writers[:xyz] = JLD2OutputWriter(model, outputs;
                                                   schedule = TimeInterval(90days),
                                                   filename = "neverworld_xyz_" * output_suffix,
                                                   with_halos = true,
                                                   overwrite_existing = true)

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(360days),
                                                        prefix = "neverworld_$(Nx)_$(Ny)_$(Nz)_checkpoint",
                                                        cleanup = true)

@info "Running..."
@show simulation

run!(simulation)

