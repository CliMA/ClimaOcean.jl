module ClimaOceanMakieExt

using Makie
using ClimaOcean.Diagnostics: simulation_report, compute_report_fields

using NumericalEarth.DataWrangling: WOAAnnual

"""
    simulation_report(ocean;
                      filename = "simulation_report.png",
                      dataset = WOAAnnual(),
                      figsize = (1800, 1200))

Generate a 3×3 diagnostic panel from the current state of `ocean` and save to `filename`.

Panels: SST, SSS, surface speed, surface vorticity, MLD,
        zonal-mean T, zonal-mean S, SST bias vs WOA, SSS bias vs WOA.
"""
function ClimaOcean.Diagnostics.simulation_report(ocean;
                                                   filename = "simulation_report.png",
                                                   dataset = WOAAnnual(),
                                                   figsize = (1800, 1200))

    f = compute_report_fields(ocean; dataset)

    fig = Figure(; size = figsize)

    # Row 1: SST, SSS, surface speed
    ax1 = Axis(fig[1, 1], title = "SST (°C)")
    hm1 = heatmap!(ax1, f.SST, colormap = :magma, colorrange = (-2, 32))
    Colorbar(fig[1, 2], hm1)

    ax2 = Axis(fig[1, 3], title = "SSS (psu)")
    hm2 = heatmap!(ax2, f.SSS, colormap = :viridis, colorrange = (32, 38))
    Colorbar(fig[1, 4], hm2)

    ax3 = Axis(fig[1, 5], title = "Surface speed (m/s)")
    hm3 = heatmap!(ax3, f.spd, colormap = :deep, colorrange = (0, 0.5))
    Colorbar(fig[1, 6], hm3)

    # Row 2: Vorticity, MLD, zonal-mean T
    ax4 = Axis(fig[2, 1], title = "Surface vorticity (s⁻¹)")
    hm4 = heatmap!(ax4, f.ζ, colormap = :balance, colorrange = (-2e-5, 2e-5))
    Colorbar(fig[2, 2], hm4)

    ax5 = Axis(fig[2, 3], title = "Mixed layer depth (m)")
    hm5 = heatmap!(ax5, f.MLD, colormap = :dense, colorrange = (0, 500))
    Colorbar(fig[2, 4], hm5)

    ax6 = Axis(fig[2, 5], title = "Zonal-mean T (°C)", xlabel = "Latitude", ylabel = "Depth (m)")
    hm6 = heatmap!(ax6, f.φ, f.z, f.T̄, colormap = :magma, colorrange = (-2, 30))
    Colorbar(fig[2, 6], hm6)

    # Row 3: Zonal-mean S, SST bias, SSS bias
    ax7 = Axis(fig[3, 1], title = "Zonal-mean S (psu)", xlabel = "Latitude", ylabel = "Depth (m)")
    hm7 = heatmap!(ax7, f.φ, f.z, f.S̄, colormap = :viridis, colorrange = (33, 37))
    Colorbar(fig[3, 2], hm7)

    ax8 = Axis(fig[3, 3], title = "SST bias vs WOA (°C)")
    hm8 = heatmap!(ax8, f.δT, colormap = :balance, colorrange = (-5, 5))
    Colorbar(fig[3, 4], hm8)

    ax9 = Axis(fig[3, 5], title = "SSS bias vs WOA (psu)")
    hm9 = heatmap!(ax9, f.δS, colormap = :balance, colorrange = (-2, 2))
    Colorbar(fig[3, 6], hm9)

    for ax in (ax1, ax2, ax3, ax4, ax5)
        hidedecorations!(ax)
    end

    save(filename, fig)
    @info "Simulation report saved to $filename"

    return fig
end

end # module
