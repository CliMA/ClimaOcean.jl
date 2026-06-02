# Figure 8: SSH RMS (m), one panel per case (Robinson, 1° lat-lon regrid).
function fig09(caches, labels, cases)
    fig = Figure(size = (800 * length(labels), 500), fontsize = 14)
    for (i, lab) in enumerate(labels)
        surface_panel!(fig, [1, 2i-1], get_field(caches[lab], :ssh_rms_latlon);
               title = "$lab: SSH RMS", colormap = :magma,
               colorrange = (0, 0.2), label = "m")
    end
    savefig(fig, "fig09_ssh_rms.png")
end
