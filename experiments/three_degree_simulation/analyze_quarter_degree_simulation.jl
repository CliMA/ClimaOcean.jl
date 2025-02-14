using Oceananigans
using ClimaOcean
using GLMakie

filename = "three_degree_simulation_surface.jld2"

# ht = FieldTimeSeries(filename, "h")
ℵt = FieldTimeSeries(filename, "ℵ")
Tt = FieldTimeSeries(filename, "T")

fig = Figure(size=(1200, 1200))

#axh = Axis(fig[1, 1])
axℵ = Axis(fig[1, 1])
axT = Axis(fig[2, 1])

Nt = length(Tt)
n = Observable(1)
#ℵn = @lift interior(ℵt[$n], :, :, 1)
ℵn = @lift ℵt[$n]
Tn = @lift interior(Tt[$n], :, :, 1)

#heatmap!(axh, hn)
heatmap!(axℵ, ℵn)
heatmap!(axT, Tn)

display(fig)

record(fig, "quarter_degree_acc.mp4", 1:Nt, framerate=12) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

