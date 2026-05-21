# Figure 1: SST - WOA bias, one panel per case (1° lat-lon regrid).
function fig01(caches, labels, cases)
    fig = Figure(size = (800 * length(labels), 500), fontsize = 14)
    for (i, lab) in enumerate(labels)
        surface_panel!(fig, [1, 2i-1], get_field(caches[lab], :sst_bias_latlon);
               title = "$lab: SST - WOA", colormap = :balance,
               colorrange = (-2.75, 2.75), label = "deg C")
    end
    savefig(fig, "fig01_sst_bias.png")
end
