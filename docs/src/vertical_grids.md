# Vertical grids

```@setup vgrids
using CairoMakie
CairoMakie.activate!(type="svg")
```

```@example vgrids
using ClimaOcean
using CairoMakie

Nz = 12
depth = 1000
zf = exponential_z_faces(; Nz, depth, scale = depth/2)

zc = [(zf[k] + zf[k+1])/2 for k in 1:Nz]
Δz = diff(zf)

fig = Figure()

axΔz1 = Axis(fig[1, 1:2]; xlabel = "z-spacing (m)", ylabel = "z (m)", title = "scale = depth /2")
axz1 = Axis(fig[1, 3])
linkaxes!(axΔz1, axz1)

scatter!(axΔz1, Δz, zc)
hidespines!(axΔz1, :t, :r)

lines!(axz1, [0, 0], [-depth, 0], color=:gray)
scatter!(axz1, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz1, 0 * zc, zc)
hidedecorations!(axz1)
hidespines!(axz1)


zf = exponential_z_faces(; Nz, depth, scale = depth/4)

zc = [(zf[k] + zf[k+1])/2 for k in 1:Nz]
Δz = diff(zf)

axΔz2 = Axis(fig[1, 4:5]; xlabel = "z-spacing (m)", ylabel = "z (m)", title = "scale = depth /4")
axz2 = Axis(fig[1, 6])
linkaxes!(axΔz2, axz2)

scatter!(axΔz2, Δz, zc)
xlims!(axΔz2, -10.2, 10.2)

# hidespines!(axΔz2, :t, :r)

lines!(axz2, [0, 0], [-depth, 0], color=:gray)
scatter!(axz2, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz2, 0 * zc, zc)
hidedecorations!(axz2)
hidespines!(axz2)

Label(fig[0, :], "z-grid")

colsize!(fig.layout, 3, Relative(0.1))
colsize!(fig.layout, 6, Relative(0.1))

fig
```
