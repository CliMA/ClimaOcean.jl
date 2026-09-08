# Figure 28: Conservation checks. Salt content is split into ocean, sea-ice, and
# ocean+ice total — the total is the quantity the freshwater mass-flux design should
# conserve; a flat total with an ocean/ice see-saw means the drift is reservoir
# exchange, a drifting total means a genuine leak (e.g. unnormalized restoring).
# Ocean heat content and ocean volume are shown alongside (volume drift = net water
# added). Water is split the same way as salt: freshwater normalization corrects only the
# atmospheric part of the surface flux, so ocean volume alone is expected to move with the
# sea-ice exchange — ocean + ice + snow, in freshwater-equivalent volume, is the quantity
# that should stay flat. Salt in Gt (1e12 kg), heat in ZJ (1e21 J), volume in km³.
function fig28(caches, labels, cases)
    fig = Figure(size = (1900, 800), fontsize = 14)

    ax_ocean  = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "Δ salt (Gt)",
                     title = "Ocean salt content drift")
    ax_ice    = Axis(fig[1, 2]; xlabel = "Time (years)", ylabel = "Δ salt (Gt)",
                     title = "Sea-ice salt content drift")
    ax_total  = Axis(fig[1, 3]; xlabel = "Time (years)", ylabel = "Δ salt (Gt)",
                     title = "Total salt content drift (ocean + ice)")
    ax_heat   = Axis(fig[2, 1]; xlabel = "Time (years)", ylabel = "Δ heat (ZJ)",
                     title = "Ocean heat content drift")
    ax_volume = Axis(fig[2, 2]; xlabel = "Time (years)", ylabel = "Δ volume (km³)",
                     title = "Ocean volume drift")
    ax_water  = Axis(fig[2, 3]; xlabel = "Time (years)", ylabel = "Δ volume (km³)",
                     title = "Total water drift (ocean + ice + snow)")

    for ax in (ax_ocean, ax_ice, ax_total, ax_volume, ax_water)
        hlines!(ax, [0]; color = (:black, 0.4), linestyle = :dash, linewidth = 1)
    end

    for (i, lab) in enumerate(labels)
        time_in_years = get_field(caches[lab], :time_in_years)
        ocean_salt    = get_field(caches[lab], :ocean_salt_content_timeseries)   ./ 1e12
        ice_salt      = get_field(caches[lab], :sea_ice_salt_content_timeseries) ./ 1e12
        total_salt    = get_field(caches[lab], :total_salt_content_timeseries)   ./ 1e12
        ocean_heat    = get_field(caches[lab], :ocean_heat_content_timeseries)   ./ 1e21
        ocean_volume  = get_field(caches[lab], :ocean_volume_timeseries)         ./ 1e9
        total_water   = get_field(caches[lab], :total_water_volume_timeseries)   ./ 1e9

        lines!(ax_ocean,  time_in_years, ocean_salt   .- ocean_salt[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
        lines!(ax_ice,    time_in_years, ice_salt     .- ice_salt[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
        lines!(ax_total,  time_in_years, total_salt   .- total_salt[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
        lines!(ax_heat,   time_in_years, ocean_heat   .- ocean_heat[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
        lines!(ax_volume, time_in_years, ocean_volume .- ocean_volume[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
        lines!(ax_water,  time_in_years, total_water  .- total_water[1];
               color = case_colors[i], linewidth = CASE_LINEWIDTH, label = lab)
    end

    Legend(fig[1:2, 4], ax_total)
    savefig(fig, "fig28_content_conservation.png")
end
