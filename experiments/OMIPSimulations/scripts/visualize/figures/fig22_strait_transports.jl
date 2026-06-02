# Figure 21: Strait transports (Bering, Drake, ITF) annual means per case.
function fig22(caches, labels, cases)
    # Bin a time series sampled at `t_seconds` into yearly means starting at year 0.
    function annual_means(t_seconds, values)
        years_full = floor.(Int, t_seconds ./ (365.25 * 86400))
        unique_years = sort(unique(years_full))
        centers = Float64[]
        means   = Float64[]
        for y in unique_years
            mask = years_full .== y
            any(mask) || continue
            push!(centers, y + 0.5)
            push!(means,   mean(values[mask]))
        end
        return centers, means
    end

    have_any = false
    for lab in labels
        st = get_field(caches[lab], :strait_transports)
        isnothing(st) || (have_any = true; break)
    end
    have_any || return

    fig = Figure(size = (1000 + 200 * length(labels), 500), fontsize = 14)
    ax_b = Axis(fig[1, 1]; xlabel = "Time (years)", ylabel = "Transport (Sv)", title = "Bering Strait")
    ax_d = Axis(fig[1, 2]; xlabel = "Time (years)", ylabel = "Transport (Sv)", title = "Drake Passage")
    ax_i = Axis(fig[1, 3]; xlabel = "Time (years)", ylabel = "Transport (Sv)", title = "Indonesian Throughflow")
    for (i, lab) in enumerate(labels)
        st = get_field(caches[lab], :strait_transports)
        isnothing(st) && continue
        tb, yb = annual_means(st.time, st.bering)
        td, yd = annual_means(st.time, st.drake)
        ti, yi = annual_means(st.time, st.itf)
        lines!(ax_b, tb, yb; color = case_colors[i], label = lab, linewidth = CASE_LINEWIDTH)
        lines!(ax_d, td, yd; color = case_colors[i], label = lab, linewidth = CASE_LINEWIDTH)
        lines!(ax_i, ti, yi; color = case_colors[i], label = lab, linewidth = CASE_LINEWIDTH)
    end
    Legend(fig[1, 4], ax_b)
    savefig(fig, "fig22_strait_transports.png")
end
