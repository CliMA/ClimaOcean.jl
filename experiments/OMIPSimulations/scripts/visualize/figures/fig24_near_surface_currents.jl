# Figure 23: Near-surface mean currents (top model cell, geographic E/N, 1° lat-lon regrid).
function fig24(caches, labels, cases)
    fig = Figure(size = (800 * length(labels), 1300), fontsize = 14)
    for (i, lab) in enumerate(labels)
        c = caches[lab]
        surface_panel!(fig, [1, 2i-1], get_field(c, :near_surface_speed_latlon);
               title = "$lab: Near-surface current speed",
               colormap = :speed, colorrange = (0, 0.5), label = "m/s")
        surface_panel!(fig, [2, 2i-1], get_field(c, :near_surface_zonal_velocity_latlon);
               title = "$lab: Zonal current (uE)",
               colormap = :balance, colorrange = (-0.5, 0.5), label = "m/s")
        surface_panel!(fig, [3, 2i-1], get_field(c, :near_surface_meridional_velocity_latlon);
               title = "$lab: Meridional current (vN)",
               colormap = :balance, colorrange = (-0.3, 0.3), label = "m/s")
    end
    savefig(fig, "fig24_near_surface_currents.png")
end
