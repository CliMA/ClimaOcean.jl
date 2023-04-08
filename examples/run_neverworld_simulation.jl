using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: min_Δx, min_Δy
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.TurbulenceClosures: HorizontalFormulation
using ClimaOcean.IdealizedSimulations: neverworld_simulation
using ClimaOcean.NearGlobalSimulations: geometric_viscosity
using Printf

minimum_turbulent_kinetic_energy = 1e-6
minimum_convective_buoyancy_flux = 1e-11
closure = CATKEVerticalDiffusivity(; minimum_turbulent_kinetic_energy,
                                   minimum_convective_buoyancy_flux)

simulation = neverworld_simulation(GPU();
                                   horizontal_size = (240, 280),
                                   longitude = (0, 60),
                                   latitude = (-70, 0),
                                   time_step = 10minutes,
                                   stop_time = 200years,
                                   closure)

model = simulation.model
grid = model.grid

@show grid

Nx, Ny, Nz = size(grid)
i = round(Int, Nx/10)
j = round(Int, Ny/2)
schedule = TimeInterval(4hours)

start_time = Ref(time_ns())

function progress(sim) 
    b = sim.model.tracers.b
    e = sim.model.tracers.e
    u, v, w = sim.model.velocities

    msg = @sprintf("Iter: %d, time: %s, extrema(b): (%6.2e, %6.2e)",
                   iteration(sim), prettytime(sim), minimum(b), maximum(b))

    msg *= @sprintf(", extrema(e): (%6.2e, %6.2e)", minimum(e), maximum(e))

    msg *= @sprintf(", max|u|: %6.2e, max|w|: %6.2e",
                    maximum(maximum(abs, q) for q in (u, v, w)), maximum(abs, w))

    try 
        κᶜ = sim.model.diffusivity_fields.Kᶜ
        msg *= @sprintf(", extrema(κᶜ): (%6.2e, %6.2e)", minimum(κᶜ), maximum(κᶜ))
    catch
    end

    elapsed = 1e-9 * (time_ns() - start_time[])
    msg *= @sprintf(", wall time: %s", prettytime(elapsed))
    start_time[] = time_ns()

    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
ζ = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)

diffusivity_fields = (; κᶜ = model.diffusivity_fields.κᶜ)

outputs = merge(model.velocities, model.tracers, diffusivity_fields)
zonally_averaged_outputs = NamedTuple(n => Average(outputs[n], dims=1) for n in keys(outputs))

simulation.output_writers[:yz] = JLD2OutputWriter(model, outputs;
                                                  schedule = TimeInterval(6days),
                                                  filename = "neverworld_yz.jld2",
                                                  indices = (i, :, :),
                                                  overwrite_existing = true)

simulation.output_writers[:zonal] = JLD2OutputWriter(model, zonally_averaged_outputs;
                                                     schedule = TimeInterval(6days),
                                                     filename = "neverworld_zonal_average.jld2",
                                                     overwrite_existing = true)

simulation.output_writers[:xy] = JLD2OutputWriter(model, outputs;
                                                  schedule = TimeInterval(6days),
                                                  filename = "neverworld_xy.jld2",
                                                  indices = (:, :, Nz),
                                                  overwrite_existing = true)

simulation.output_writers[:xyz] = JLD2OutputWriter(model, outputs;
                                                   schedule = TimeInterval(90days),
                                                   filename = "neverworld_xyz.jld2",
                                                   overwrite_existing = true)

run!(simulation)

#=
using GLMakie

bxyt = FieldTimeSeries("neverworld_xy.jld2", "b")
uxyt = FieldTimeSeries("neverworld_xy.jld2", "u")
ζxyt = FieldTimeSeries("neverworld_xy.jld2", "ζ")
exyt = FieldTimeSeries("neverworld_xy.jld2", "e")
κxyt = FieldTimeSeries("neverworld_xy.jld2", "κᶜ")

byzt = FieldTimeSeries("neverworld_yz.jld2", "b")
uyzt = FieldTimeSeries("neverworld_yz.jld2", "u")
αyzt = FieldTimeSeries("neverworld_yz.jld2", "α")
ζyzt = FieldTimeSeries("neverworld_yz.jld2", "ζ")
eyzt = FieldTimeSeries("neverworld_yz.jld2", "e")
κyzt = FieldTimeSeries("neverworld_yz.jld2", "κᶜ")

t = bxyt.times
Nt = length(t)

fig = Figure(resolution=(1800, 600))

axbxy = Axis(fig[1, 1], title="b(x, y)")
axuxy = Axis(fig[1, 2], title="u(x, y)")
axζxy = Axis(fig[1, 3], title="ζ(x, y)")
axexy = Axis(fig[1, 4], title="e(x, y)")
axκxy = Axis(fig[1, 5], title="κᶜ(x, y)")

axbyz = Axis(fig[2, 1], title="b(y, z)")
axuyz = Axis(fig[2, 2], title="| ∂z u |²(y, z)")
axζyz = Axis(fig[2, 3], title="ζ(y, z)")
axeyz = Axis(fig[2, 4], title="e(y, z)")
axκyz = Axis(fig[2, 5], title="κᶜ(y, z)")

slider = Slider(fig[2, 1:5], range=1:Nt, startvalue=Nt)
n = slider.value

bxyn = @lift interior(bxyt[$n], :, :, 1)
uxyn = @lift interior(uxyt[$n], :, :, 1)
ζxyn = @lift interior(ζxyt[$n], :, :, 1)
exyn = @lift interior(exyt[$n], :, :, 1)
κxyn = @lift interior(κxyt[$n], :, :, 1)

byzn = @lift interior(byzt[$n], 1, :, :)
uyzn = @lift interior(αyzt[$n], 1, :, :)
ζyzn = @lift interior(ζyzt[$n], 1, :, :)
eyzn = @lift interior(eyzt[$n], 1, :, :)
κyzn = @lift interior(κyzt[$n], 1, :, :)

ulim = 1.0
ζlim = 5e-5
κlim = 100
x, y, z = nodes(bxyt)

heatmap!(axbxy, bxyn, colorrange=(-0.06, 0.06))
heatmap!(axuxy, uxyn, colormap=:balance, colorrange=(-ulim, ulim))
heatmap!(axζxy, ζxyn, colormap=:balance, colorrange=(-ζlim, ζlim))
heatmap!(axexy, exyn, colormap=:solar, colorrange=(1e-6, 1e-2))
heatmap!(axκxy, κxyn, colormap=:solar, colorrange=(0, κlim))

heatmap!(axbyz, y, z, byzn, colorrange=(0.0, 0.06))
heatmap!(axuyz, y, z, uyzn, colormap=:balance, colorrange=(-ulim, ulim))
heatmap!(axζyz, y, z, ζyzn, colormap=:balance, colorrange=(-ζlim, ζlim))
heatmap!(axeyz, y, z, eyzn, colormap=:solar, colorrange=(1e-6, 1e-2))
heatmap!(axκyz, y, z, κyzn, colormap=:solar, colorrange=(0, κlim))

# display(fig)

record(fig, "catke_neverworld.mp4", 1:Nt, framerate=48) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end
=#
