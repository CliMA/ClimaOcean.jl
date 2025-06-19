# Vertical grids

A few vertical grids are implemented within the [VerticalGrids](@ref ClimaOcean.VerticalGrids) module.

### Exponential spacing

This grid implements a vertical grid with an exponentially growing spacing.
The faces are distributed in the range ``[-L_z, 0]`` using the scaling

```math
\frac{\exp{[(z + L_z) / h]} - 1}{\exp{(L_z / h)} - 1}
```

which for ``h / L \to \infty`` it reduces to

```math
1 + z / L_z
```

The above implies that the vertical spacings grow linearly with depth, with the smallest spacing being near the surface.
For ``h / L \to \infty``, the grid becomes uniform.

To showcase how the scale ``h`` affects the grid, we construct below two such exponential grids,
one with ``h / L = 1/5`` and the second with ``h / L = 1/2``.


```@setup vgrids
using CairoMakie
CairoMakie.activate!(type="svg")
```

```@example vgrids
using ClimaOcean

Nz = 10
depth = 1000

scale = depth / 5

zf  = exponential_z_faces(; Nz, depth, scale) # z-faces
zc = [(zf[k] + zf[k+1])/2 for k in 1:Nz]     # z-centers
Δz = diff(zf)                                # spacing between z-faces

using CairoMakie

fig = Figure()

axΔz1 = Axis(fig[1, 1]; xlabel = "z-spacing (m)", ylabel = "z (m)", title = "scale = depth / 5")
axz1 = Axis(fig[1, 2])
linkaxes!(axΔz1, axz1)

scatter!(axΔz1, Δz, zc)
hidespines!(axΔz1, :t, :r)

lines!(axz1, [0, 0], [-depth, 0], color=:gray)
scatter!(axz1, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz1, 0 * zc, zc)
hidedecorations!(axz1)
hidespines!(axz1)


scale = depth / 2
zf = exponential_z_faces(; Nz, depth, scale) # z-faces
zc = [(zf[k] + zf[k+1])/2 for k in 1:Nz]     # z-centers
Δz = diff(zf)                                # spacing between z-faces

axΔz2 = Axis(fig[1, 3]; xlabel = "z-spacing (m)", ylabel = "z (m)", title = "scale = depth / 2")
axz2 = Axis(fig[1, 4])
linkaxes!(axΔz2, axz2)

scatter!(axΔz2, Δz, zc)
hidespines!(axΔz2, :t, :r)

lines!(axz2, [0, 0], [-depth, 0], color=:gray)
scatter!(axz2, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz2, 0 * zc, zc)
hidedecorations!(axz2)
hidespines!(axz2)


Label(fig[0, :], "z-grid", fontsize = 18)

colsize!(fig.layout, 2, Relative(0.1))
colsize!(fig.layout, 4, Relative(0.1))

fig
```

Both grid spacings grow linearly and when summed up give the total depth.
But with the larger ``h / L_z`` is, the smaller the rate the spacings increase with depth becomes.

A ridiculously large value of ``h / L_z`` (approaching infinity) gives a uniform grid:

```@example vgrids
exponential_z_faces(; Nz, depth, scale = 1e20*depth)
```
