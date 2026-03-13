# # [Sixth-degree tripolar ocean–sea ice simulation](@id sixth-degree-distributed)
#
# This page visualizes the output from a 1/6° coupled ocean–sea ice simulation
# on a `TripolarGrid` (2160×1080), distributed across 4 GPUs and forced by
# repeat-year JRA55 atmospheric reanalysis for two years.
#
# The simulation is run separately via MPI:
# ```
# srun -n 4 julia --project examples/sixth_degree_tripolar_ocean_sea_ice.jl
# ```
# and the output JLD2 files are loaded here for visualization.

using CairoMakie
using Oceananigans
using Oceananigans.Units
using Printf

# ### Load saved output

To = FieldTimeSeries("sixth_degree_coupled_ocean_surface.jld2", "T"; backend = OnDisk())
uo = FieldTimeSeries("sixth_degree_coupled_ocean_surface.jld2", "u"; backend = OnDisk())
vo = FieldTimeSeries("sixth_degree_coupled_ocean_surface.jld2", "v"; backend = OnDisk())

hi = FieldTimeSeries("sixth_degree_coupled_sea_ice.jld2", "h"; backend = OnDisk())
ℵi = FieldTimeSeries("sixth_degree_coupled_sea_ice.jld2", "ℵ"; backend = OnDisk())

times = To.times
Nt = length(times)
n = Observable(Nt)

# ### Land mask

land = interior(To.grid.immersed_boundary.bottom_height) .≥ 0

# ### Observables for plotting

Toₙ = @lift begin
    Tₙ = interior(To[$n])
    Tₙ[land] .= NaN
    view(Tₙ, :, :, 1)
end

# Compute surface speed from u and v.
uoₙ = Field{Face, Center, Nothing}(uo.grid)
voₙ = Field{Center, Face, Nothing}(vo.grid)
so = Field(sqrt(uoₙ^2 + voₙ^2))

soₙ = @lift begin
    parent(uoₙ) .= parent(uo[$n])
    parent(voₙ) .= parent(vo[$n])
    compute!(so)
    sₙ = interior(so)
    sₙ[land] .= NaN
    view(sₙ, :, :, 1)
end

# Effective sea ice thickness (h × ℵ).
heₙ = @lift begin
    hₙ = interior(hi[$n])
    ℵₙ = interior(ℵi[$n])
    hₙ[land] .= NaN
    view(hₙ, :, :, 1) .* view(ℵₙ, :, :, 1)
end

# ### Snapshot

fig = Figure(size = (1200, 1000))

title = @lift string("1/6° distributed simulation after ", prettytime(times[$n] - times[1]))

axT = Axis(fig[1, 1])
axs = Axis(fig[1, 3])
axh = Axis(fig[2, 1])

hmT = heatmap!(axT, Toₙ, colorrange = (-1, 32), colormap = :magma,  nan_color = :lightgray)
hms = heatmap!(axs, soₙ, colorrange = (0, 0.5),  colormap = :deep,   nan_color = :lightgray)
hmh = heatmap!(axh, heₙ, colorrange = (0, 4),     colormap = :blues,  nan_color = :lightgray)

Colorbar(fig[1, 2], hmT, label = "Surface Temperature (°C)")
Colorbar(fig[1, 4], hms, label = "Surface Speed (m s⁻¹)")
Colorbar(fig[2, 2], hmh, label = "Effective ice thickness (m)")

for ax in (axT, axs, axh)
    hidedecorations!(ax)
end

Label(fig[0, :], title)

save("sixth_degree_snapshot.png", fig)
nothing #hide

# ![](sixth_degree_snapshot.png)

# ### Movie

CairoMakie.record(fig, "sixth_degree_simulation.mp4", 1:Nt, framerate = 8) do nn
    n[] = nn
end
nothing #hide

# ![](sixth_degree_simulation.mp4)
