# Figure 6: Sea-ice concentration bias (model − HadISST1, 1979–2007) on
# polar stereographic projections (45°N–90°N and 45°S–90°S) for March and
# September. Mirrors Adcroft et al. (2019) Figure 25: diverging colormap
# centred at zero, with the observed 15% SIC contour overlaid in black
# on every panel.
# Rows: NH March, SH March, NH September, SH September. Columns: cases.
function fig06(caches, labels, cases)
    ncases = length(labels)
    fig = Figure(size = (500 * ncases, 500 * 4), fontsize = 14)

    function maybe_plot(row_idx, lab_idx, lab, bias_sym, obs_sym, hemisphere, season)
        bias = get_field(caches[lab], bias_sym)
        obs  = get_field(caches[lab], obs_sym)
        if isnothing(bias) || isnothing(obs)
            @warn "fig06: skipping $lab $season ($hemisphere) — bias or obs is nothing. Did the HadISST download succeed?"
            return
        end
        polar_panel!(fig, [row_idx, 2*lab_idx - 1], bias;
                     hemisphere,
                     title = "$lab: $season SIC bias ($(hemisphere == :north ? "NH" : "SH"))",
                     colormap = Reverse(:RdBu), colorrange = (-1, 1),
                     label = "model − HadISST",
                     obs_contour = obs,
                     obs_levels = [0.15], obs_color = :black, obs_linewidth = 2.5)
    end

    for (i, lab) in enumerate(labels)
        maybe_plot(1, i, lab, :sic_march_bias_latlon,     :sic_march_obs_latlon,     :north, "March")
        maybe_plot(2, i, lab, :sic_march_bias_latlon,     :sic_march_obs_latlon,     :south, "March")
        maybe_plot(3, i, lab, :sic_september_bias_latlon, :sic_september_obs_latlon, :north, "September")
        maybe_plot(4, i, lab, :sic_september_bias_latlon, :sic_september_obs_latlon, :south, "September")
    end
    savefig(fig, "fig06_seaice_conc_bias.png")
end
