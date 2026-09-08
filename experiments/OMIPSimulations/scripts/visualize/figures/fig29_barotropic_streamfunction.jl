# Figure 29: Barotropic streamfunction (1° lat-lon regrid), globally and over the North Atlantic.
# Ψ = -∫∫ u dz dy with Ψ = 0 at the southern boundary, so subtropical gyres are positive, subpolar
# gyres negative, and the ACC negative with magnitude equal to the Drake Passage transport. The North
# Atlantic panel is the diagnostic for the subpolar gyre and the NAC path; observational estimates put
# the subpolar gyre near -30 Sv and the North Atlantic subtropical gyre near +30 Sv.
function fig29(caches, labels, cases)
    fig = Figure(size = (800 * length(labels), 900), fontsize = 14)
    for (i, lab) in enumerate(labels)
        Ψ = get_field(caches[lab], :barotropic_streamfunction_latlon)

        surface_panel!(fig, [1, 2i-1], Ψ;
               title = "$lab: Barotropic streamfunction",
               colormap = :balance, colorrange = (-60, 60), label = "Sv")

        north_atlantic_panel!(fig, [2, 2i-1], Ψ;
               title = "$lab: North Atlantic",
               colormap = :balance, colorrange = (-40, 40), label = "Sv")
    end
    savefig(fig, "fig29_barotropic_streamfunction.png")
end

# Subpolar/subtropical North Atlantic window. Subsets the data before projecting, the same way
# `polar_panel!` clips its cap, so the region fills the panel instead of sitting in a global view.
function north_atlantic_panel!(fig, pos, data; kwargs...)
    lonlims = (-80.0, 10.0)
    latlims = (35.0, 70.0)
    i_keep = findall(λ -> lonlims[1] <= λ <= lonlims[2], to_minus180_180(LATLON_LON_CENTERS, data)[1])
    j_keep = findall(φ -> latlims[1] <= φ <= latlims[2], LATLON_LAT_CENTERS)
    lon_shifted, data_shifted = to_minus180_180(LATLON_LON_CENTERS, data)
    return geo_panel!(fig, pos, data_shifted[i_keep, j_keep];
                      x = lon_shifted[i_keep],
                      y = LATLON_LAT_CENTERS[j_keep],
                      projection = "+proj=stere +lat_0=55 +lon_0=-35 +lat_ts=55",
                      lonlims, latlims,
                      kwargs...)
end
