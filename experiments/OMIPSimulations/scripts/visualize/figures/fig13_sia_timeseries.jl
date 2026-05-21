# Figure 12: Sea-ice area time series (Arctic + Antarctic) per case.
function fig13(caches, labels, cases)
    m2_to_million_km2 = 1e-12
    fig = Figure(size = (600 + 200 * length(labels), 500), fontsize = 14)
    ax_arctic = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "SIA (Million km²)",
                     title = "Arctic sea-ice area")
    for (i, lab) in enumerate(labels)
        ice = get_field(caches[lab], :sea_ice_diagnostics)
        time_in_years = [Dates.value(d - ice.snapshot_dates[1]) / (365.25 * 86400 * 1000)
                         for d in ice.snapshot_dates]
        lines!(ax_arctic, time_in_years, ice.arctic_area .* m2_to_million_km2;
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    ax_antarctic = Axis(fig[1, 2]; xlabel = "Time (years)", ylabel = "SIA (Million km²)",
                        title = "Antarctic sea-ice area")
    for (i, lab) in enumerate(labels)
        ice = get_field(caches[lab], :sea_ice_diagnostics)
        time_in_years = [Dates.value(d - ice.snapshot_dates[1]) / (365.25 * 86400 * 1000)
                         for d in ice.snapshot_dates]
        lines!(ax_antarctic, time_in_years, ice.antarctic_area .* m2_to_million_km2;
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    Legend(fig[1, 3], ax_arctic)
    savefig(fig, "fig13_sia_timeseries.png")
end
