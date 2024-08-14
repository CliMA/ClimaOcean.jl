using ClimaOcean
using CairoMakie
using Oceananigans
using Oceananigans.Units
using Printf

time_indices = Colon()
Qswt = ClimaOcean.JRA55.JRA55_field_time_series(:downwelling_shortwave_radiation; time_indices)
rht = ClimaOcean.JRA55.JRA55_field_time_series(:relative_humidity; time_indices)

function lonlat2xyz(lons::AbstractVector, lats::AbstractVector)
    x = [cosd(lat) * cosd(lon) for lon in lons, lat in lats]
    y = [cosd(lat) * sind(lon) for lon in lons, lat in lats]
    z = [sind(lat)             for lon in lons, lat in lats]
    return (x, y, z)
end

function lonlat2xyz(lons::AbstractMatrix, lats::AbstractMatrix)
    x = @. cosd(lats) * cosd(lons)
    y = @. cosd(lats) * sind(lons)
    z = @. sind(lats)
    return (x, y, z)
end

λ, φ, z = nodes(Qswt)
x, y, z = lonlat2xyz(λ, φ)

times = Qswt.times
Nt = length(times)
n = Observable(1)

Qswn = @lift interior(Qswt[$n], :, :, 1)
rhn = @lift interior(rht[$n], :, :, 1)

fig = Figure(size=(1400, 700))

axsw = Axis3(fig[1, 1], aspect=(1, 1, 1))
axrh = Axis3(fig[1, 2], aspect=(1, 1, 1))

label = @lift string("Repeat-year JRA55 forcing on year-day ",
                     @sprintf("%.1f", times[$n] / days))

Label(fig[0, 1:2], label, fontsize=24)

sf = surface!(axsw, x, y, z, color=Qswn, colorrange=(0, 1200))
Colorbar(fig[2, 1], sf,
         vertical = false,
         width = Relative(0.5),
         flipaxis = false,
         label = "Downwelling shortwave radiation (W m⁻²)")

sf = surface!(axrh, x, y, z, color=rhn, colormap=:grays, colorrange=(0, 100))
Colorbar(fig[2, 2], sf,
         vertical = false,
         width = Relative(0.5),
         flipaxis = false,
         label = "Relative humidity (%)")

colgap!(fig.layout, 1, Relative(-0.15))
rowgap!(fig.layout, 1, Relative(-0.2))
rowgap!(fig.layout, 2, Relative(-0.2))

for ax in (axsw, axrh)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.viewmode = :fit # so that the sphere does not zoom in and out while rotating
end

display(fig)

snapshot_interval = 3hours  # JRE55 time-resolution
rotation_period = 60days
rotation_rate = 2π / rotation_period

record(fig, "JRA55_data.mp4", 1:Nt, framerate=16) do nn
    @info nn/Nt

    n[] = nn

    for ax in (axsw, axrh)
        ax.azimuth = nn * snapshot_interval * rotation_rate
    end
end
