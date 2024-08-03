# # A quick look at ECCO data
#
# ClimaOcean can download and utilize data from the "ECCO" state estimate,
# which stands for "Estimating the Circulation and Climate of the Ocean" --- two!
#
# This script shows how to download three-dimensional temperature and salinity fields
# from ECCO, and makes a short animation to showcase the fields' content.
#
# For this example we need Oceananigans for Field utilities, CairoMakie for plotting
# Printf for nice labeling, and of course ClimaOcean to actually download and construct
# the ECCO fields.

using Oceananigans
using CairoMakie
using Printf

using ClimaOcean: ECCO

# The function `ecco_field` provided by `ClimaOcean.DataWrangling.ECCO` will automatically
# download ECCO data if it doesn't already exist at the default location.

T = ECCO.ecco_field(:temperature)
S = ECCO.ecco_field(:salinity)

# Next, we massage the ECCO data by inserting NaNs in "land cells", which
# are diagnosed by having an unphysically low temperature.

Tp = parent(T)
Sp = parent(S)
Sp[Tp .< -10] .= NaN
Tp[Tp .< -10] .= NaN

# # Plotting ECCO data
#
# We're ready to plot. We'll make an animation
# that depicts how the ECCO data changes with depth.

fig = Figure(size=(1200, 1400))

axT = Axis(fig[1, 1])
axS = Axis(fig[2, 1])

λ, φ, z = nodes(T)

# To make an animation that scrolls through the 3D temperature
# and salinity fields, we make an Observable for the vertical index,
# and then construct slices of T, S using the Observable, `k`.

grid = T.grid
Nz = size(grid, 3)
k = Observable(Nz)

Tk = @lift interior(T, :, :, $k)
Sk = @lift interior(S, :, :, $k)

# Finally, we make a nice plot with a label that displays depth, colorbars,
# and light gray within land cells.

hmT = heatmap!(axT, λ, φ, Tk, nan_color=:lightgray, colorrange=(-2, 30), colormap=:thermal)
hmS = heatmap!(axS, λ, φ, Sk, nan_color=:lightgray, colorrange=(31, 37), colormap=:haline)

Colorbar(fig[1, 2], hmT, label="Temperature (ᵒC)")
Colorbar(fig[2, 2], hmS, label="Salinity (psu)")

z = znodes(grid, Center())
depth_str = @lift @sprintf("%.1f meters depth", -z[$k])
text!(axT, 50, 50, text=depth_str, justification=:center, fontsize=24)
text!(axS, 50, 50, text=depth_str, justification=:center, fontsize=24)

# # Making the animation
#
# This animation is a little fancy. We start by displaying the surface
# field, then we scroll through depth to the bottom and pause there.
# Next, we scroll back to the surface and pause.

stillframes = 10
movingframes = Nz

record(fig, "ECCO_temperature_salinity.mp4", framerate=4) do io

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
nothing #hide

# ![](ECCO_temperature_salinity.mp4)
