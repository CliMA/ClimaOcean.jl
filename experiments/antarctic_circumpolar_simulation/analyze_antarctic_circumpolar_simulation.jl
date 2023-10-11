using Oceananigans
using GLMakie

filename = "antarctic_circumpolar_surface.jld2"

ut = FieldTimeSeries(filename, "u")
Tt = FieldTimeSeries(filename, "T")
St = FieldTimeSeries(filename, "S")

times = Tt.times
Nt = length(times)

function set_zero_to_NaN!(u)
    up = parent(u)
    up[up .== 0] .= NaN
    return nothing
end

for n = 1:Nt
    un = ut[n]
    Tn = Tt[n]
    Sn = St[n]
    set_zero_to_NaN!(Tn)
    set_zero_to_NaN!(Sn)
    if n > 1
        set_zero_to_NaN!(un)
    end
end

set_theme!(Theme(fontsize=32))
fig = Figure(resolution=(2600, 2000))

axT = Axis(fig[1, 1], xlabel="Longitude (degrees)", ylabel="Latitude (degrees)")
axS = Axis(fig[2, 1], xlabel="Longitude (degrees)", ylabel="Latitude (degrees)")
axu = Axis(fig[3, 1], xlabel="Longitude (degrees)", ylabel="Latitude (degrees)")

slider = Slider(fig[4, 1], range=1:Nt, startvalue=1)
n = slider.value

Tn = @lift interior(Tt[$n], :, :, 1)
Sn = @lift interior(St[$n], :, :, 1)
un = @lift interior(ut[$n], :, :, 1)

λc, φc, zc = nodes(Tt)
λu, φc, zc = nodes(ut)

hm = heatmap!(axT, λc, φc, Tn, colorrange=(-1, 30), nan_color=:gray, colormap=:thermal)
Colorbar(fig[1, 2], hm, label="Temperature (ᵒC)")

hm = heatmap!(axS, λc, φc, Sn, colorrange=(31, 35), nan_color=:gray, colormap=:haline)
Colorbar(fig[2, 2], hm, label="Salinity (psu)")

hm = heatmap!(axu, λu, φc, un, colorrange=(-0.2, 0.2),  nan_color=:gray, colormap=:redblue)
Colorbar(fig[3, 2], hm, label="Zonal velocity (m s⁻¹)")



display(fig)

