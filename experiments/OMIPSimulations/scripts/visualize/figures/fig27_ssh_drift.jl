# Figure 27: Global-mean free-surface (zosga) time series — the Boussinesq
# mass-conservation drift that shifts the absolute SSH. Plotted absolute (not
# demeaned) so a nonzero global mean is visible at a glance.
function fig27(caches, labels, cases)
    fig = Figure(size = (600 + 200 * length(labels), 450), fontsize = 14)
    ax = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "⟨η⟩ (m)",
              title = "Global-mean free-surface displacement")
    hlines!(ax, [0]; color = (:black, 0.4), linestyle = :dash, linewidth = 1)
    for (i, lab) in enumerate(labels)
        ssh           = get_field(caches[lab], :global_mean_ssh_timeseries)
        time_in_years = get_field(caches[lab], :time_in_years)
        lines!(ax, time_in_years, ssh;
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    Legend(fig[1, 2], ax)
    savefig(fig, "fig27_ssh_drift.png")
end
