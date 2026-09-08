# Figure 30: Arctic freshwater export through Fram and Davis Straits, liquid and solid, annual means.
#
# The test of the hosing hypothesis: if the Labrador/Irminger fresh bias arrives from the north, the
# export here is too strong. Sign follows `strait_transports` — positive northward, so export out of
# the Arctic is negative and a too-strong export plots *below* the observational band.
#
# Reference values are southward exports from Serreze et al. (2006), "The large-scale freshwater
# cycle of the Arctic" (JGR 111, C11010), and for Davis Strait from Curry et al. (2014), JPO 44,
# 1244-1266, both relative to a 34.8 reference salinity — the same one the diagnostic uses.
const ARCTIC_FRESHWATER_OBSERVATIONS = (
    fram  = (liquid = (-2660.0, 500.0), solid = (-2300.0, 500.0)),
    davis = (liquid = (-2930.0, 400.0), solid = ( -315.0, 150.0)),
)

function fig30(caches, labels, cases)
    # Bin a time series sampled at `t_seconds` into yearly means starting at year 0.
    function annual_means(t_seconds, values)
        years_full = floor.(Int, t_seconds ./ (365.25 * 86400))
        centers = Float64[]
        means   = Float64[]
        for y in sort(unique(years_full))
            mask = years_full .== y
            any(mask) || continue
            push!(centers, y + 0.5)
            push!(means,   mean(values[mask]))
        end
        return centers, means
    end

    have_any = false
    for lab in labels
        fw = get_field(caches[lab], :strait_freshwater_transports)
        isnothing(fw) || (have_any = true; break)
    end
    have_any || return

    fig = Figure(size = (1400, 900), fontsize = 14)
    axes = Dict{Tuple{Symbol, Symbol}, Any}()

    for (row, strait) in enumerate((:fram, :davis)), (col, phase) in enumerate((:liquid, :solid))
        title = string(uppercasefirst(string(strait)), " Strait — ", string(phase))
        ax = Axis(fig[row, col]; xlabel = "Time (years)",
                  ylabel = "Freshwater transport (km³ yr⁻¹)", title)
        axes[(strait, phase)] = ax

        # Observational band. Negative is southward export, so a model that over-exports sits below.
        centre, spread = getproperty(getproperty(ARCTIC_FRESHWATER_OBSERVATIONS, strait), phase)
        hlines!(ax, [centre]; color = OBS_COLOR, linestyle = OBS_LINESTYLE, linewidth = OBS_LINEWIDTH,
                label = "Observed")
        hspan!(ax, centre - spread, centre + spread; color = (OBS_COLOR, 0.15))
        hlines!(ax, [0.0]; color = :black, linewidth = 0.8)
    end

    for (i, lab) in enumerate(labels)
        fw = get_field(caches[lab], :strait_freshwater_transports)
        isnothing(fw) && continue
        for strait in (:fram, :davis), phase in (:liquid, :solid)
            component = getproperty(fw, phase)
            haskey(component, strait) || continue
            t, y = annual_means(component.time, getproperty(component, strait))
            lines!(axes[(strait, phase)], t, y;
                   color = case_colors[i], label = lab, linewidth = CASE_LINEWIDTH)
        end
    end

    Legend(fig[1:2, 3], axes[(:fram, :liquid)], merge = true, unique = true)
    savefig(fig, "fig30_arctic_freshwater.png")
end
