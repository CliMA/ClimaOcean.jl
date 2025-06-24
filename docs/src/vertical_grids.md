# Vertical grids

A few vertical grids are implemented within the [VerticalGrids](@ref ClimaOcean.VerticalGrids) module.

### Exponential spacing

The [`ExponentialFaces`](@ref) method returns a vertical grid with faces at depths following an exponential profile.
The faces are distributed in the range ``[-L, 0]`` using the scaling

```math
\frac{\exp{[(z + L) / h]} - 1}{\exp{(L / h)} - 1}
```

which varies from 0 (at ``z = -L``) to 1 (at the surface, ``z = 0``).

The exponential scaling above implies that the vertical spacings grow linearly with depth at a rate inversely proportional to ``h / L``, with the smallest spacing being near the surface.

At the limit ``h / L \to \infty`` the scaling reduces to ``1 + z/L`` and thus the grid becomes uniformly spaced.


```@setup vgrids
using CairoMakie
set_theme!(Theme(Lines = (linewidth = 3,)))
CairoMakie.activate!(type="svg")
```

```@example vgrids
using ClimaOcean
using ClimaOcean.VerticalGrids: exponential_profile

using CairoMakie

depth = 1000
z = range(-depth, stop=0, length=501)

fig = Figure()
ax = Axis(fig[1, 1], ylabel="z/L")

scale = depth / 20
lines!(ax, exponential_profile.(z, depth, scale), z / depth, label="h / L = $(scale / depth)")
scale = depth / 5
lines!(ax, exponential_profile.(z, depth, scale), z / depth, label="h / L = $(scale / depth)")
scale = depth / 2
lines!(ax, exponential_profile.(z, depth, scale), z / depth, label="h / L = $(scale / depth)")
scale = 1e20 * depth
lines!(ax, exponential_profile.(z, depth, scale), z / depth, label="h / L → ∞")

axislegend(ax, position=:rb)

fig
```

Let's see how [`ExponentialFaces`](@ref) works:

```@example vgrids
using ClimaOcean

Nz = 10
depth = 1000

z = ExponentialFaces(Nz, depth; scale)
```

The above returns a callable object which gives that ``k``-th face, e.g.,

```@example vgrids
[z(k) for k in 1:Nz+1]
```

To showcase how the scale ``h`` affects the grid, we construct below two such exponential grids,
one with ``h / L = 1/5`` and the second with ``h / L = 1/2``.

```@example vgrids
using ClimaOcean
using Oceananigans

Nz = 10
depth = 1000

scale = depth / 5
z = ExponentialFaces(Nz, depth; scale)
grid = RectilinearGrid(; size=Nz, z, topology=(Flat, Flat, Bounded))
zf = znodes(grid, Face())
zc = znodes(grid, Center())
Δz = zspacings(grid, Center())
Δz = view(Δz, 1, 1, :)  # for plotting

using CairoMakie

fig = Figure()

axΔz1 = Axis(fig[1, 1]; xlabel = "z-spacing (m)", ylabel = "z (m)", title = "scale = depth / 5")
axz1 = Axis(fig[1, 2])
linkaxes!(axΔz1, axz1)

lΔz = lines!(axΔz1, - zf * (depth/scale) / Nz, zf, color=(:grey, 0.5))
scatter!(axΔz1, Δz, zc)
hidespines!(axΔz1, :t, :r)

lines!(axz1, [0, 0], [-depth, 0], color=:gray)
scatter!(axz1, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz1, 0 * zc, zc)
hidedecorations!(axz1)
hidespines!(axz1)


scale = depth / 2
z = ExponentialFaces(Nz, depth; scale)
grid = RectilinearGrid(; size=Nz, z, topology=(Flat, Flat, Bounded))
zf = znodes(grid, Face())
zc = znodes(grid, Center())
Δz = zspacings(grid, Center())
Δz = view(Δz, 1, 1, :)  # for plotting

axΔz2 = Axis(fig[1, 3]; xlabel = "z-spacing (m)", ylabel = "z (m)", title = "scale = depth / 2")
axz2 = Axis(fig[1, 4])
linkaxes!(axΔz2, axz2)

lΔz = lines!(axΔz2, - zf * (depth/scale) / Nz, zf, color=(:grey, 0.5))
scatter!(axΔz2, Δz, zc)
hidespines!(axΔz2, :t, :r)

lines!(axz2, [0, 0], [-depth, 0], color=:gray)
scatter!(axz2, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz2, 0 * zc, zc)
hidedecorations!(axz2)
hidespines!(axz2)


colsize!(fig.layout, 2, Relative(0.1))
colsize!(fig.layout, 4, Relative(0.1))

legend = Legend(fig[2, :], [lΔz], ["slope = (depth / scale) / Nz"], orientation = :horizontal)

fig
```

For both grids, the spacings grow linearly with depth and sum up to the total depth.
But with the larger ``h / L`` is, the smaller the rate is that the spacings increase with depth.

A ridiculously large value of ``h / L`` (approximating infinity) gives a uniform grid:

```@example vgrids
z = ExponentialFaces(Nz, depth, scale = 1e20*depth)
[z(k) for k in 1:Nz+1]
```

A downside of the above ais that we don't have tight control on the minimum spacing at the surface.
To prescribe the surface spacing we need to play around with the scale ``h`` and the number of vertical cells ``N_z``.

### Stretched ``z`` faces

The [`stretched_vertical_faces`](@ref) method allows a tighter control on the vertical spacing at the surface.
That is, we can prescribe a constant spacing over the top `surface_layer_height`  below which the grid spacing
increases following a prescribed stretching law.
The downside is that neither the final grid depth nor the total number of vertical cells can be prescribed.
The final depth we get is greater or equal from what we prescribed via the keyword argument `depth`
and also the total number of faces depends on the stretching law.

The three grids below have constant 20-meter spacing for the top 120 meters.
We prescribe to all `depth = 750` meters and we apply power-law stretching for depths below 120 meters.
The bigger the power-law stretching factor is, the further the last face goes beyond the prescribed depth and/or with less total number of cells.

```@example vgrids
depth = 750
surface_layer_Δz = 20
surface_layer_height = 120

z = StretchedFaces(; depth,
                   surface_layer_Δz,
                   surface_layer_height,
                   stretching = PowerLawStretching(1.08))
grid = RectilinearGrid(; size=length(z), z, topology=(Flat, Flat, Bounded))
zf = znodes(grid, Face())
zc = znodes(grid, Center())
Δz = zspacings(grid, Center())
Δz = view(Δz, 1, 1, :)  # for plotting

fig = Figure(size=(800, 550))

axΔz1 = Axis(fig[1, 1];
             xlabel = "z-spacing (m)",
             ylabel = "z (m)",
             title = "PowerLawStretching(1.09)\n $(length(zf)-1) cells\n bottom face at z = $(zf[1]) m\n ")

axz1 = Axis(fig[1, 2])

ldepth = hlines!(axΔz1, -depth, color = :salmon, linestyle=:dash)
lzbottom = hlines!(axΔz1, zf[1], color = :grey)
scatter!(axΔz1, Δz, zc)
hidespines!(axΔz1, :t, :r)

lines!(axz1, [0, 0], [zf[1], 0], color=:gray)
scatter!(axz1, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz1, 0 * zc, zc)
hidedecorations!(axz1)
hidespines!(axz1)


z = StretchedFaces(; depth,
                   surface_layer_Δz,
                   surface_layer_height,
                   stretching = PowerLawStretching(1.04))
grid = RectilinearGrid(; size=length(z), z, topology=(Flat, Flat, Bounded))
zf = znodes(grid, Face())
zc = znodes(grid, Center())
Δz = zspacings(grid, Center())
Δz = view(Δz, 1, 1, :)  # for plotting

axΔz2 = Axis(fig[1, 3];
             xlabel = "z-spacing (m)",
             ylabel = "z (m)",
             title = "PowerLawStretching(1.04)\n $(length(zf)-1) cells\n bottom face at z = $(zf[1]) m\n ")
axz2 = Axis(fig[1, 4])

ldepth = hlines!(axΔz2, -depth, color = :salmon, linestyle=:dash)
lzbottom = hlines!(axΔz2, zf[1], color = :grey)
scatter!(axΔz2, Δz, zc)
hidespines!(axΔz2, :t, :r)

lines!(axz2, [0, 0], [zf[1], 0], color=:gray)
scatter!(axz2, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz2, 0 * zc, zc)
hidedecorations!(axz2)
hidespines!(axz2)


z = StretchedFaces(; depth,
                   surface_layer_Δz,
                   surface_layer_height,
                   stretching = PowerLawStretching(1.04),
                   constant_bottom_spacing_depth = 500)
grid = RectilinearGrid(; size=length(z), z, topology=(Flat, Flat, Bounded))
zf = znodes(grid, Face())
zc = znodes(grid, Center())
Δz = zspacings(grid, Center())
Δz = view(Δz, 1, 1, :)  # for plotting

axΔz3 = Axis(fig[1, 5];
             xlabel = "z-spacing (m)",
             ylabel = "z (m)",
             title = "PowerLawStretching(1.04)\n $(length(zf)-1) cells\n bottom face at z = $(zf[1]) m\n constant spacing below 500 m")
axz3 = Axis(fig[1, 6])

ldepth = hlines!(axΔz3, -depth, color = :salmon, linestyle=:dash)
lzbottom = hlines!(axΔz3, zf[1], color = :grey)
scatter!(axΔz3, Δz, zc)

hidespines!(axΔz3, :t, :r)

lines!(axz3, [0, 0], [zf[1], 0], color=:gray)
scatter!(axz3, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz3, 0 * zc, zc)
hidedecorations!(axz3)
hidespines!(axz3)


linkaxes!(axΔz1, axz1, axΔz2, axz2, axΔz3, axz3)

Legend(fig[2, :], [ldepth, lzbottom], ["prescribed depth", "bottom-most z-face"], orientation = :horizontal)

colsize!(fig.layout, 2, Relative(0.1))
colsize!(fig.layout, 4, Relative(0.1))
colsize!(fig.layout, 6, Relative(0.1))

fig
```
