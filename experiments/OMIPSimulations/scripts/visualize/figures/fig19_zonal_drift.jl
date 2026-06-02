# Figure 18: Zonal-mean drift (T - WOA, S - WOA, b - b₀) as filled contours.
function fig19(caches, labels, cases)
    latitude = zonal_latitude_centers()
    temperature_bias_levels = range(-3, 3; length = 13)
    salinity_bias_levels    = range(-0.75, 0.75; length = 13)

    # Buoyancy drift levels: symmetric around 0; magnitude set from data so
    # the colormap actually resolves the field across cases.
    buoyancy_drift_max = maximum(lab -> begin
                       v = filter(isfinite, vec(get_field(caches[lab], :zonal_buoyancy_drift)))
                       isempty(v) ? 0.0 : maximum(abs, v)
                   end, labels)
    buoyancy_drift_max = buoyancy_drift_max == 0 ? 1e-3 : buoyancy_drift_max
    buoyancy_drift_levels = range(-buoyancy_drift_max, buoyancy_drift_max; length = 13)

    fig = Figure(size = (600 * length(labels), 900), fontsize = 14)
    for (i, lab) in enumerate(labels)
        c = caches[lab]
        depth = get_field(c, :depth)

        ax = Axis(fig[1, 2i-1]; xlabel = "Latitude", ylabel = "Depth (m)",
                  title = "$lab: Zonal T - WOA")
        hm = contourf!(ax, latitude, depth, get_field(c, :zonal_temperature_bias);
                       levels = temperature_bias_levels, colormap = :balance,
                       extendlow = :auto, extendhigh = :auto)
        Colorbar(fig[1, 2i], hm; label = "deg C"); ylims!(ax, (-5500, 0))

        ax = Axis(fig[2, 2i-1]; xlabel = "Latitude", ylabel = "Depth (m)",
                  title = "$lab: Zonal S - WOA")
        hm = contourf!(ax, latitude, depth, get_field(c, :zonal_salinity_bias);
                       levels = salinity_bias_levels, colormap = :balance,
                       extendlow = :auto, extendhigh = :auto)
        Colorbar(fig[2, 2i], hm; label = "PSU"); ylims!(ax, (-5500, 0))

        ax = Axis(fig[3, 2i-1]; xlabel = "Latitude", ylabel = "Depth (m)",
                  title = "$lab: Zonal b - b(t=0)")
        hm = contourf!(ax, latitude, depth, get_field(c, :zonal_buoyancy_drift);
                       levels = buoyancy_drift_levels, colormap = :balance,
                       extendlow = :auto, extendhigh = :auto)
        Colorbar(fig[3, 2i], hm; label = "m/s²"); ylims!(ax, (-5500, 0))
    end
    savefig(fig, "fig19_zonal_drift.png")
end
