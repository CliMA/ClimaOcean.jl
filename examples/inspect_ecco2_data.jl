using Oceananigans
using ClimaOcean: ECCO2
using CairoMakie
using Printf

T = ECCO2.ecco2_field(:temperature)
S = ECCO2.ecco2_field(:salinity)

grid = T.grid
Nz = size(grid, 3)

# Put NaNs at unphysical (probably, under land) temperature values
Tp = parent(T)
Sp = parent(S)
Sp[Tp .< -10] .= NaN
Tp[Tp .< -10] .= NaN

fig = Figure(resolution=(1200, 1400))

axT = Axis(fig[1, 1])
axS = Axis(fig[2, 1])

λ, φ, z = nodes(T)

k = Observable(Nz)
Tk = @lift interior(T, :, :, $k)
Sk = @lift interior(S, :, :, $k)

hmT = heatmap!(axT, λ, φ, Tk, nan_color=:lightgray, colorrange=(-2, 30), colormap=:thermal)
hmS = heatmap!(axS, λ, φ, Sk, nan_color=:lightgray, colorrange=(31, 37), colormap=:haline)

Colorbar(fig[1, 2], hmT, label="Temperature (ᵒC)")
Colorbar(fig[2, 2], hmS, label="Salinity (psu)")

z = znodes(grid, Center())
depth_str = @lift @sprintf("%.1f meters depth", -z[$k])
text!(axT, 50, 50, text=depth_str, justification=:center, fontsize=24)
text!(axS, 50, 50, text=depth_str, justification=:center, fontsize=24)

stillframes = 10
movingframes = Nz

record(fig, "ECCO2_temperature_salinity.gif", framerate=4) do io

    [recordframe!(io) for _ = 1:stillframes]

    for kk in Nz:-4:1
        k[] = kk
        recordframe!(io)
    end

    [recordframe!(io) for _ = 1:stillframes]

    for kk in 1:4:Nz
        k[] = kk
        recordframe!(io)
    end

    [recordframe!(io) for _ = 1:stillframes]
end

