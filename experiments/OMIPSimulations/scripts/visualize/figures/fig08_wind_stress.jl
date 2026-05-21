# Figure 7: Wind stress (zonal, meridional) and (optional) NCEP biases (1° lat-lon regrid).
function fig08(caches, labels, cases)
    has_ncep = any(lab -> !isnothing(get_field(caches[lab], :zonal_wind_stress_bias_ncep)), labels)
    nrows = has_ncep ? 4 : 2
    fig = Figure(size = (800 * length(labels), 450 * nrows), fontsize = 14)
    for (i, lab) in enumerate(labels)
        c = caches[lab]
        surface_panel!(fig, [1, 2i-1], get_field(c, :zonal_wind_stress_latlon);
               title = "$lab: Zonal wind stress", colormap = :balance,
               colorrange = (-0.3, 0.3), label = "N/m²")
        surface_panel!(fig, [2, 2i-1], get_field(c, :meridional_wind_stress_latlon);
               title = "$lab: Meridional wind stress", colormap = :balance,
               colorrange = (-0.3, 0.3), label = "N/m²")
        if has_ncep
            zonal_bias      = get_field(c, :zonal_wind_stress_bias_ncep_latlon)
            meridional_bias = get_field(c, :meridional_wind_stress_bias_ncep_latlon)
            !isnothing(zonal_bias) && surface_panel!(fig, [3, 2i-1], zonal_bias;
                title = "$lab: τx - NCEP", colormap = :balance,
                colorrange = (-0.15, 0.15), label = "N/m²")
            !isnothing(meridional_bias) && surface_panel!(fig, [4, 2i-1], meridional_bias;
                title = "$lab: τy - NCEP", colormap = :balance,
                colorrange = (-0.15, 0.15), label = "N/m²")
        end
    end
    savefig(fig, "fig08_wind_stress.png")
end
