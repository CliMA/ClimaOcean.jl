
using GLMakie

fig = Figure(resolution = (800, 600))
axQ = Axis(fig[1, 1], title = "Heat Flux")
axS = Axis(fig[1, 2], title = "Salinity Flux")
axu = Axis(fig[2, 1], title = "Zonal stress")
axv = Axis(fig[2, 2], title = "Meridional stress")

ρ  = reference_density(earth.model.ocean.model)
cp = heat_capacity(earth.model.ocean.model)

heatmap!(axQ, stats.Jᵀ.avg * ρ * cp, colormap = :balance)
heatmap!(axS, stats.Jˢ.avg, colormap = :balance)
heatmap!(axu, stats.τx.avg, colormap = :balance)
heatmap!(axv, stats.τy.avg, colormap = :balance)