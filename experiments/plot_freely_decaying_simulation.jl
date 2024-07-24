using CairoMakie
using Oceananigans
using JLD2

sea_ice_filename = "freely_decaying_regional_simulation_heat_only_sea_ice_thickness.jld2"
fields_filename =  "freely_decaying_regional_simulation_heat_only_fields.jld2"

ht = FieldTimeSeries(sea_ice_filename, "h")
Tt = FieldTimeSeries(fields_filename, "T")
times = ht.times
Nt = length(times)

for n = 1:Nt
    Tn = interior(Tt[n])
    hn = interior(ht[n])
    land = Tn .== 0
    Tn[land] .= NaN
    hn[land] .= NaN
end

fig = Figure(resolution=(1200, 600))
axT = Axis(fig[1, 1], xlabel="Longitude", ylabel="Latitude")
axh = Axis(fig[2, 1], xlabel="Longitude", ylabel="Latitude")

n = Observable(1)
Tn = @lift interior(Tt[$n], :, :, 1)
hn = @lift interior(ht[$n], :, :, 1)

λ, φ, z = nodes(Tt)

hmT = heatmap!(axT, λ, φ, Tn, nan_color=:lightyellow, colormap=:thermal, colorrange=(-2, 22))
hmh = heatmap!(axh, λ, φ, hn, nan_color=:lightyellow, colormap=:grays, colorrange=(0, 1))

Colorbar(fig[1, 2], hmT, label="Temperature (ᵒC)")
Colorbar(fig[2, 2], hmh, label="Sea ice thickness (m)")

record(fig, "free_decay_southern_ocean.mp4", 1:Nt, framerate=12) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end
