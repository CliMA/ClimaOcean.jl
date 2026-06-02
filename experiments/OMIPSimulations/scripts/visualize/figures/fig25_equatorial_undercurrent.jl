# Figure 24: Equatorial undercurrent sections, after Ringler et al. (2013) fig. 5.
#   * Meridional section of zonal velocity at 140°W (lat ∈ [-8°, 10°])
#   * Zonal section of zonal velocity along the equator (lon ∈ [143°E, 95°W])
# Observations (Johnson et al. 2002) are stubbed until the file is wired up.

function fig25(caches, labels, cases)
    ncases  = length(labels)
    levels  = collect(-100:10:100)
    heavy   = collect(-100:50:100)
    has_obs = any(lab -> !isnothing(get_field(caches[lab], :euc_obs_meridional)) ||
                         !isnothing(get_field(caches[lab], :euc_obs_equatorial)), labels)
    nrows   = has_obs ? 4 : 2

    fig = Figure(size = (700 * ncases, 500 * nrows), fontsize = 14)

    function draw_meridional!(row, col, title, lats, depth, data)
        ax = Axis(fig[row, 2col - 1]; xlabel = "Latitude", ylabel = "Depth (m)", title)
        hm = heatmap!(ax, lats, depth, data;
                      colormap = :balance, colorrange = (-100, 100), nan_color = :lightgray)
        contour!(ax, lats, depth, data;
                 levels, color = :black, linewidth = 0.5)
        contour!(ax, lats, depth, data;
                 levels = heavy, color = :black, linewidth = 1.4)
        xlims!(ax, (-8, 10)); ylims!(ax, (-400, 0))
        Colorbar(fig[row, 2col], hm; label = "cm/s")
    end

    function draw_equatorial!(row, col, title, lons, depth, data)
        ax = Axis(fig[row, 2col - 1]; xlabel = "Longitude (°E)", ylabel = "Depth (m)", title)
        hm = heatmap!(ax, lons, depth, data;
                      colormap = :balance, colorrange = (-100, 100), nan_color = :lightgray)
        contour!(ax, lons, depth, data;
                 levels, color = :black, linewidth = 0.5)
        contour!(ax, lons, depth, data;
                 levels = heavy, color = :black, linewidth = 1.4)
        xlims!(ax, (143, 265)); ylims!(ax, (-400, 0))
        Colorbar(fig[row, 2col], hm; label = "cm/s")
    end

    for (i, lab) in enumerate(labels)
        c     = caches[lab]
        depth = get_field(c, :euc_depth)

        # Model sections (m/s → cm/s)
        uE_m_model   = 100 .* get_field(c, :euc_meridional_section)
        lats_m       = get_field(c, :euc_meridional_lats)
        uE_eq_model  = 100 .* get_field(c, :euc_equatorial_section)
        lons_eq      = get_field(c, :euc_equatorial_lons)

        model_row_m  = has_obs ? 3 : 1
        model_row_eq = has_obs ? 4 : 2

        # Observation panels (top half) — only if any case provides them.
        if has_obs
            obs_m  = get_field(c, :euc_obs_meridional)
            obs_eq = get_field(c, :euc_obs_equatorial)
            if !isnothing(obs_m)
                draw_meridional!(1, i, "$lab: u(140°W) Johnson 2002",
                                 obs_m.lats, obs_m.depth, obs_m.uE)
            end
            if !isnothing(obs_eq)
                draw_equatorial!(2, i, "$lab: u(equator) Johnson 2002",
                                 obs_eq.lons, obs_eq.depth, obs_eq.uE)
            end
        end

        draw_meridional!(model_row_m,  i, "$lab: u(140°W) model",  lats_m,  depth, uE_m_model)
        draw_equatorial!(model_row_eq, i, "$lab: u(equator) model", lons_eq, depth, uE_eq_model)
    end

    savefig(fig, "fig25_equatorial_undercurrent.png")
end
