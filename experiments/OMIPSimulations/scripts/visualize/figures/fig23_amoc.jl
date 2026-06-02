# Figure 22: Atlantic Meridional Overturning Circulation streamfunction.
# ψ_atl(j, z) computed per j-row on the model grid (no regridding) from
# the saved `vvol = v·Aʸ` field, summed zonally over the Atlantic basin
# mask and cumulatively integrated from the bottom up. Sv = 1e6 m³/s.
function fig23(caches, labels, cases)
    levels = -5:17
    fig = Figure(size = (600 * length(labels), 500), fontsize = 14)
    for (i, lab) in enumerate(labels)
        c     = caches[lab]
        depth = get_field(c, :depth)
        lats  = get_field(c, :amoc_latitudes)
        ψ     = get_field(c, :amoc)
        ψ   .+= 0.0
    
                 valid = 2:length(lats)-1
                    lats  = lats[valid]
                      ψ     = ψ[valid, :]
        
        @show size(depth), size(lats), size(ψ)
        ax = Axis(fig[1, 2i-1]; xlabel = "Latitude", ylabel = "Depth (m)",
                  title = "$lab: AMOC ψ")
        hm = contourf!(ax, lats, depth, ψ;
                       levels, colormap = :balance,
                       extendlow = :auto, extendhigh = :auto)
        contour!(ax, lats, depth, ψ; levels, color = :black, linewidth = 0.6)
        Colorbar(fig[1, 2i], hm; label = "Sv")
        ylims!(ax, (-5500, 0))
        xlims!(ax, (-35, 65))
    end
    savefig(fig, "fig23_amoc.png")
end
