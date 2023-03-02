using ClimaOcean.NearGlobalSimulations: aqua_planet_simulation
using Oceananigans
using Oceananigans.Units
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
using GLMakie

# closure = RiBasedVerticalDiffusivity() 
longitude = (-30, 30)
simulation = aqua_planet_simulation(GPU(), size=(240, 480, 30))

model = simulation.model
grid = model.grid
Nz = grid.Nz
slices_save_interval = 1day
output_prefix = "aqua_planet"
dir = "."

slices = Dict(:surface => (:, :, Nz))

for (name, indices) in slices
    outputs = Dict()

    for name in keys(model.tracers)
        c = model.tracers[name]
        outputs[name] = Field(c; indices)
    end

    outputs[:u] = Field(model.velocities.u; indices)
    outputs[:v] = Field(model.velocities.v; indices)
    outputs[:w] = Field(model.velocities.w; indices)
    outputs[:η] = model.free_surface.η
    outputs[:ζ] = VerticalVorticityField(model.grid, model.velocities; indices)

    simulation.output_writers[name] = JLD2OutputWriter(model, outputs; dir,
                                                       schedule = TimeInterval(slices_save_interval),
                                                       filename = output_prefix * "_$name",
                                                       with_halos = true,
                                                       overwrite_existing = true)
end


@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

simulation.Δt = 1minutes
simulation.stop_time = 300days
run!(simulation)

@info "Simulation took $(prettytime(simulation.run_wall_time))."

ut = FieldTimeSeries("aqua_planet_surface.jld2", "u")
vt = FieldTimeSeries("aqua_planet_surface.jld2", "v")
Tt = FieldTimeSeries("aqua_planet_surface.jld2", "T")
et = FieldTimeSeries("aqua_planet_surface.jld2", "e")

times = ut.times
Nt = length(times)

fig = Figure(resolution=(1600, 1000))

slider = Slider(fig[1, 1:2], range=1:Nt, startvalue=1)
n = slider.value

un = @lift interior(ut[$n], :, :, 1)
vn = @lift interior(vt[$n], :, :, 1)
Tn = @lift interior(Tt[$n], :, :, 1)
en = @lift interior(et[$n], :, :, 1)

ax = Axis(fig[2, 1])
heatmap!(ax, un)

ax = Axis(fig[2, 2])
heatmap!(ax, vn)

ax = Axis(fig[3, 1])
heatmap!(ax, Tn)

ax = Axis(fig[3, 2])
heatmap!(ax, en)

display(fig)

