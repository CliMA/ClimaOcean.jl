# Figure 3: Time-mean SSH (row 1) and SSH - ECCO bias (row 2), 1° lat-lon regrid.
function fig03(caches, labels, cases)
    fig = Figure(size = (800 * length(labels), 900), fontsize = 14)
    for (i, lab) in enumerate(labels)
        surface_panel!(fig, [1, 2i-1], get_field(caches[lab], :ssh_latlon);
               title = "$lab: Time-mean SSH", colormap = :balance,
               colorrange = (-1.5, 1.5), label = "m")
        surface_panel!(fig, [2, 2i-1], get_field(caches[lab], :ssh_bias_ecco_latlon);
               title = "$lab: SSH - ECCO (1992–2012), demeaned",
               colormap = :balance, colorrange = (-0.5, 0.5), label = "m")
    end
    savefig(fig, "fig03_ssh.png")
end
