# Figure 13: Arctic sea-ice volume time series per case.
function fig14(caches, labels, cases)
    m3_to_thousand_km3 = 1e-12
    fig = Figure(size = (400 + 200 * length(labels), 500), fontsize = 14)
    ax = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "Ice volume (10³ km³)",
              title = "Arctic sea-ice volume")
    for (i, lab) in enumerate(labels)
        ice = get_field(caches[lab], :sea_ice_diagnostics)
        time_in_years = [Dates.value(d - ice.snapshot_dates[1]) / (365.25 * 86400 * 1000)
                         for d in ice.snapshot_dates]
        lines!(ax, time_in_years, ice.arctic_volume .* m3_to_thousand_km3;
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    Legend(fig[1, 2], ax)
    savefig(fig, "fig14_arctic_volume_timeseries.png")
end
