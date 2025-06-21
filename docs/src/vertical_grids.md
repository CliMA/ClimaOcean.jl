# Vertical grids

A few vertical grids are implemented within the [VerticalGrids](@ref ClimaOcean.VerticalGrids) module.

### Exponential spacing

The [`exponential_vertical_faces`](@ref) method returns a vertical grid with an exponentially growing spacing.
The faces are distributed in the range ``[-L_z, 0]`` using the scaling

```math
\frac{\exp{[(z + L_z) / h]} - 1}{\exp{(L_z / h)} - 1}
```

which varies from 0 (at ``z = -L_z``) to 1 (at the surface, ``z = 0``).
At the limit ``h / L \to \infty`` the scaling above reduces to

```math
1 + z / L_z
```

The above implies that the vertical spacings grow linearly with depth at a rate inversely proportional to ``h / L_z ``, with the smallest spacing being near the surface.
For ``h / L \to \infty``, the grid thus becomes uniform.

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
zgrid = exponential_vertical_faces(; Nz, depth, scale)

zf = z_faces(zgrid)       # z-faces
zc = z_centers(zgrid)     # z-centers
Δz = diff(zf)             # spacing between z-faces

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
zgrid = exponential_vertical_faces(; Nz, depth, scale)

zf = z_faces(zgrid)       # z-faces
zc = z_centers(zgrid)     # z-centers
Δz = diff(zf)             # spacing between z-faces

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


Label(fig[0, :], "z-grid", fontsize = 18)

colsize!(fig.layout, 2, Relative(0.1))
colsize!(fig.layout, 4, Relative(0.1))

legend = Legend(fig[2, :], [lΔz], ["slope = (depth / scale) / Nz"], orientation = :horizontal)

fig
```

Both grid spacings grow linearly and when summed up give the total depth.
But with the larger ``h / L_z`` is, the smaller the rate the spacings increase with depth becomes.

A ridiculously large value of ``h / L_z`` (approaching infinity) gives a uniform grid:

```@example vgrids
z_faces(exponential_vertical_faces(; Nz, depth, scale = 1e20*depth))
```

A downside of the above grid is that we don't have tight control on the minimum spacing at the surface.
To prescribe the surface spacing we need to play around with the scale ``h`` and the number of grid points ``N_z``.

### Stretched ``z`` faces

The [`stretched_vertical_faces`](@ref) method allows a tighter control on the vertical spacing at the surface.
That is, we can have a constant spacing over the top `surface_layer_height` and then the grid spacing
increase following a prescribed stretching law.
The downside is that neither the final grid depth nor the number of faces can be prescribed.
The final depth we get is greater or equal from what we prescribed via the keyword argument `depth`
and also the total number of faces depends on the stretching law.

The two grids below have constant 20-meter spacing for the top 120 meters.
We prescribe to both `depth = 750` meters and we apply power-law stretching for depths below 120 meters.
The bigger the power-law stretching factor is, the further the last face goes beyond the prescribed depth and/or with less total number of cells.

```@example vgrids
depth = 750
surface_layer_Δz = 20
surface_layer_height = 120

zgrid = stretched_vertical_faces(; depth,
                                 surface_layer_Δz,
                                 surface_layer_height,
                                 stretching = PowerLawStretching(1.15))
zf = z_faces(zgrid)       # z-faces
zc = z_centers(zgrid)     # z-centers
Δz = diff(zf)             # spacing between z-faces

fig = Figure()

axΔz1 = Axis(fig[1, 1];
             xlabel = "z-spacing (m)",
             ylabel = "z (m)",
             title = "PowerLawStretching(1.15)\n $(length(zf)-1) cells\n bottom face at z = $(zf[1]) m")

axz1 = Axis(fig[1, 2])

scatter!(axΔz1, Δz, zc)
hidespines!(axΔz1, :t, :r)

lines!(axz1, [0, 0], [zf[1], 0], color=:gray)
scatter!(axz1, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz1, 0 * zc, zc)
hidedecorations!(axz1)
hidespines!(axz1)


zgrid = stretched_vertical_faces(; depth,
                                 surface_layer_Δz,
                                 surface_layer_height,
                                 stretching = PowerLawStretching(1.05))
zf = z_faces(zgrid)       # z-faces
zc = z_centers(zgrid)     # z-centers
Δz = diff(zf)             # spacing between z-faces

axΔz2 = Axis(fig[1, 3];
             xlabel = "z-spacing (m)",
             ylabel = "z (m)",
             title = "PowerLawStretching(1.05)\n $(length(zf)-1) cells\n bottom face at z = $(zf[1]) m")
axz2 = Axis(fig[1, 4])

scatter!(axΔz2, Δz, zc)
hidespines!(axΔz2, :t, :r)

lines!(axz2, [0, 0], [zf[1], 0], color=:gray)
scatter!(axz2, 0 * zf, zf, marker=:hline, color=:gray, markersize=20)
scatter!(axz2, 0 * zc, zc)
hidedecorations!(axz2)
hidespines!(axz2)

linkaxes!(axΔz1, axz1, axΔz2, axz2)
Label(fig[0, :], "z-grid", fontsize = 18)

colsize!(fig.layout, 2, Relative(0.1))
colsize!(fig.layout, 4, Relative(0.1))

fig
```
