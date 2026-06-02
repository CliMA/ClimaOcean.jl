# Figure 2: SSS - WOA bias, one panel per case (1° lat-lon regrid).
function fig02(caches, labels, cases)
    fig = Figure(size = (800 * length(labels), 500), fontsize = 14)
    for (i, lab) in enumerate(labels)
        surface_panel!(fig, [1, 2i-1], get_field(caches[lab], :sss_bias_latlon);
               title = "$lab: SSS - WOA", colormap = :balance,
               colorrange = (-3, 3), label = "PSU")
    end
    savefig(fig, "fig02_sss_bias.png")
end
