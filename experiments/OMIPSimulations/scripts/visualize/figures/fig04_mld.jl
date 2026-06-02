# Figure 4: MLD seasonal min/max with optional dBM reference row (1° lat-lon regrid).
function fig04(caches, labels, cases)
    ncases = length(labels)
    label_with_dbm = findfirst(lab -> !isnothing(get_field(caches[lab], :mld_min_dbm)), labels)
    nrows = if isnothing(label_with_dbm)
        2
    elseif ncases >= 2
        3
    else
        4
    end
    fig = Figure(size = (800 * ncases, 450 * nrows), fontsize = 14)
    for (i, lab) in enumerate(labels)
        surface_panel!(fig, [1, 2i-1], get_field(caches[lab], :mld_min_latlon);
               title = "$lab: Min MLD (summer)",
               colormap = Reverse(:deep), colorrange = (0, 70), label = "m")
        surface_panel!(fig, [2, 2i-1], get_field(caches[lab], :mld_max_latlon);
               title = "$lab: Max MLD (winter)",
               colormap = Reverse(:deep), colorrange = (0, 500), label = "m")
    end
    if !isnothing(label_with_dbm)
        ref_label = labels[label_with_dbm]
        min_pos = [3, 1]
        max_pos = ncases >= 2 ? [3, 3] : [4, 1]
        surface_panel!(fig, min_pos, get_field(caches[ref_label], :mld_min_dbm_latlon);
               title = "dBM climatology: Min MLD",
               colormap = Reverse(:deep), colorrange = (0, 70), label = "m")
        surface_panel!(fig, max_pos, get_field(caches[ref_label], :mld_max_dbm_latlon);
               title = "dBM climatology: Max MLD",
               colormap = Reverse(:deep), colorrange = (0, 500), label = "m")
    end
    savefig(fig, "fig04_mld.png")
end
