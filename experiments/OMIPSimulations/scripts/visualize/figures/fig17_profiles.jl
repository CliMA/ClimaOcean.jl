# Figure 16: Time-mean horizontal-mean T and S profiles.
function fig17(caches, labels, cases)
    fig = Figure(size = (600 + 200 * length(labels), 600), fontsize = 14)
    ax_temperature = Axis(fig[1, 1]; xlabel = "Temperature (deg C)", ylabel = "Depth (m)",
                          title = "Horizontal-mean temperature")
    for (i, lab) in enumerate(labels)
        lines!(ax_temperature,
               get_field(caches[lab], :horizontal_mean_temperature_profile),
               get_field(caches[lab], :depth);
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    ylims!(ax_temperature, (-5500, 0))
    ax_salinity = Axis(fig[1, 2]; xlabel = "Salinity (PSU)", ylabel = "Depth (m)",
                       title = "Horizontal-mean salinity")
    for (i, lab) in enumerate(labels)
        lines!(ax_salinity,
               get_field(caches[lab], :horizontal_mean_salinity_profile),
               get_field(caches[lab], :depth);
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    ylims!(ax_salinity, (-5500, 0))
    Legend(fig[1, 3], ax_temperature)
    savefig(fig, "fig17_profiles.png")
end
