# Figure 6: Net heat flux (W/m², row 1) and net freshwater flux (kg/m²/s, row 2)
# (1° lat-lon regrid).
function fig07(caches, labels, cases)
    fig = Figure(size = (800 * length(labels), 900), fontsize = 14)
    for (i, lab) in enumerate(labels)
        surface_panel!(fig, [1, 2i-1], get_field(caches[lab], :heat_flux_latlon);
               title = "$lab: Net heat flux", colormap = :balance,
               colorrange = (-200, 200), label = "W/m^2")
        surface_panel!(fig, [2, 2i-1], get_field(caches[lab], :freshwater_flux_latlon);
               title = "$lab: Net freshwater flux", colormap = :balance,
               colorrange = (-1e-5, 1e-5), label = "kg/m^2/s")
    end
    savefig(fig, "fig07_surface_fluxes.png")
end
