# Figure 15: Global-mean T and S drift time series.
function fig16(caches, labels, cases)
    fig = Figure(size = (600 + 200 * length(labels), 450), fontsize = 14)
    ax_temperature = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "ΔT (deg C)",
                          title = "Global-mean temperature drift")
    for (i, lab) in enumerate(labels)
        temperature   = get_field(caches[lab], :global_mean_temperature_timeseries)
        time_in_years = get_field(caches[lab], :time_in_years)
        lines!(ax_temperature, time_in_years, temperature .- temperature[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    ax_salinity = Axis(fig[1, 2]; xlabel = "Time (years)", ylabel = "ΔS (PSU)",
                       title = "Global-mean salinity drift")
    for (i, lab) in enumerate(labels)
        salinity      = get_field(caches[lab], :global_mean_salinity_timeseries)
        time_in_years = get_field(caches[lab], :time_in_years)
        lines!(ax_salinity, time_in_years, salinity .- salinity[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    Legend(fig[1, 3], ax_temperature)
    savefig(fig, "fig16_drift.png")
end
