# Figure 11: Arctic sea-ice volume climatology (model + PIOMAS).
function fig12(caches, labels, cases)
    month_names = ["J","F","M","A","M","J","J","A","S","O","N","D"]
    m3_to_thousand_km3 = 1e-12
    obs = piomas_monthly()

    fig = Figure(size = (400 + 200 * length(labels), 500), fontsize = 14)
    ax = Axis(fig[1, 1]; xlabel = "Month", ylabel = "Ice volume (10³ km³)",
              title = "Arctic sea-ice volume", xticks = (1:12, month_names))
    if !isnothing(obs)
        band!(ax, 1:12,
              obs.volume_monthly .- obs.volume_monthly_std,
              obs.volume_monthly .+ obs.volume_monthly_std;
              color = (OBS_COLOR, 0.25))
        lines!(ax, 1:12, obs.volume_monthly;
            color = OBS_COLOR, linewidth = OBS_LINEWIDTH, linestyle = OBS_LINESTYLE, label = "PIOMAS")
    end
    for (i, lab) in enumerate(labels)
        lines!(ax, 1:12,
               get_field(caches[lab], :sea_ice_diagnostics).arctic_volume_monthly .* m3_to_thousand_km3;
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end
    Legend(fig[1, 2], ax)
    savefig(fig, "fig12_arctic_volume.png")
end
