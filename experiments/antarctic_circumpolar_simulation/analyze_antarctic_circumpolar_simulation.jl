using Oceananigans
using GLMakie

filename = "antarctic_circumpolar_surface.jld2"

ut = FieldTimeSeries(filename, "u")
Tt = FieldTimeSeries(filename, "T")
St = FieldTimeSeries(filename, "S")

times = Tt.times
Nt = length(times)

fig = Figure()
axT = Axis(fig[1, 1])
axS = Axis(fig[2, 1])
axu = Axis(fig[3, 1])

slider = Slider(fig[4, 1], range=1:Nt, startvalue=1)
n = slider.value

Tn = @lift interior(Tt[$n], :, :, 1)
Sn = @lift interior(St[$n], :, :, 1)
un = @lift interior(ut[$n], :, :, 1)

heatmap!(axT, Tn)
heatmap!(axS, Sn)
heatmap!(axu, un)

display(fig)

