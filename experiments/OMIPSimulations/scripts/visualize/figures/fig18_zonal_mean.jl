# Figure 17: Zonal-mean T, S, b sections per case (with WOA / initial overlays).
function fig18(caches, labels, cases)
    latitude = zonal_latitude_centers()
    temperature_levels = -2:2:30
    salinity_levels    = 33:0.25:37
    buoyancy_levels    = range(-0.04, 0.02, length = 13)

    fig = Figure(size = (600 * length(labels), 1200), fontsize = 14)
    for (i, lab) in enumerate(labels)
        c = caches[lab]
        depth = get_field(c, :depth)

        ax = Axis(fig[1, 2i-1]; xlabel = "Latitude", ylabel = "Depth (m)",
                  title = "$lab: Zonal T")
        hm = heatmap!(ax, latitude, depth, get_field(c, :zonal_temperature);
                      colormap = :thermal, colorrange = (-2, 30), nan_color = :lightgray)
        contour!(ax, latitude, depth, get_field(c, :zonal_woa_temperature);
                 levels = temperature_levels, color = :white, linestyle = :dash, linewidth = 1.2)
        contour!(ax, latitude, depth, get_field(c, :zonal_temperature);
                 levels = temperature_levels, color = :black, linewidth = 0.8)
        Colorbar(fig[1, 2i], hm; label = "deg C"); ylims!(ax, (-5500, 0))

        ax = Axis(fig[2, 2i-1]; xlabel = "Latitude", ylabel = "Depth (m)",
                  title = "$lab: Zonal S")
        hm = heatmap!(ax, latitude, depth, get_field(c, :zonal_salinity);
                      colormap = :haline, colorrange = (33, 37), nan_color = :lightgray)
        contour!(ax, latitude, depth, get_field(c, :zonal_woa_salinity);
                 levels = salinity_levels, color = :white, linestyle = :dash, linewidth = 1.2)
        contour!(ax, latitude, depth, get_field(c, :zonal_salinity);
                 levels = salinity_levels, color = :black, linewidth = 0.8)
        Colorbar(fig[2, 2i], hm; label = "PSU"); ylims!(ax, (-5500, 0))

        ax = Axis(fig[3, 2i-1]; xlabel = "Latitude", ylabel = "Depth (m)",
                  title = "$lab: Zonal b")
        hm = heatmap!(ax, latitude, depth, get_field(c, :zonal_buoyancy);
                      colormap = :balance, nan_color = :lightgray)
        contour!(ax, latitude, depth, get_field(c, :zonal_initial_buoyancy);
                 levels = buoyancy_levels, color = :white, linestyle = :dash, linewidth = 1.2)
        contour!(ax, latitude, depth, get_field(c, :zonal_buoyancy);
                 levels = buoyancy_levels, color = :black, linewidth = 0.8)
        Colorbar(fig[3, 2i], hm; label = "m/s²"); ylims!(ax, (-5500, 0))
    end
    savefig(fig, "fig18_zonal_mean.png")
end
