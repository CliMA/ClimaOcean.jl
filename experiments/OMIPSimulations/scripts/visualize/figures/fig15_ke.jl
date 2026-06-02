# Figure 14: Global-mean kinetic energy time series.
function fig15(caches, labels, cases)
    fig = Figure(size = (600 + 150 * length(labels), 450), fontsize = 14)
    ax = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "KE (m²/s²)",
              title = "Global-mean kinetic energy")
    for (i, lab) in enumerate(labels)
        kinetic_energy = get_field(caches[lab], :kinetic_energy)
        isempty(kinetic_energy) && continue
        time_in_years  = get_field(caches[lab], :kinetic_energy_time_in_years)
        lines!(ax, time_in_years, kinetic_energy;
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    Legend(fig[1, 2], ax)
    savefig(fig, "fig15_ke.png")
end
