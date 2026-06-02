# Figure 5: Sea-ice concentration, split into NH (polar stereographic, 45°N–90°N)
# and SH (polar stereographic, 45°S–90°S) for both March and September.
# Rows: NH March, SH March, NH September, SH September. Columns: cases.
function fig05(caches, labels, cases)
    ncases = length(labels)
    fig = Figure(size = (500 * ncases, 500 * 4), fontsize = 14)

    function maybe_plot(row_idx, lab_idx, lab, sym, hemisphere, season)
        data = get_field(caches[lab], sym)
        isnothing(data) && return
        polar_panel!(fig, [row_idx, 2*lab_idx - 1], data;
                     hemisphere,
                     title = "$lab: $season ($(hemisphere == :north ? "NH" : "SH"))",
                     colormap = :ice, colorrange = (0, 1), label = "fraction")
    end

    for (i, lab) in enumerate(labels)
        maybe_plot(1, i, lab, :sic_march_latlon,     :north, "March")
        maybe_plot(2, i, lab, :sic_march_latlon,     :south, "March")
        maybe_plot(3, i, lab, :sic_september_latlon, :north, "September")
        maybe_plot(4, i, lab, :sic_september_latlon, :south, "September")
    end
    savefig(fig, "fig05_seaice_conc.png")
end
