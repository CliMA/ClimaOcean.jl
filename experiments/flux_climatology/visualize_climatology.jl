
using GLMakie

mask = ClimaOcean.ECCO.ECCO_mask()
mask = Array(interior(mask, :, :, 50)) |> Array{Float64}
mask[mask .== 1] .= NaN

fig = Figure(resolution = (800, 600))
axQ = Axis(fig[1, 1], title = "Heat Flux")
axS = Axis(fig[1, 2], title = "Salinity Flux")
axu = Axis(fig[2, 1], title = "Zonal stress")
axv = Axis(fig[2, 2], title = "Meridional stress")

ρ  = reference_density(earth.model.ocean.model)
cp = heat_capacity(earth.model.ocean.model)

Qavg = compute!(Field(stats.Jᵀ.avg * ρ * cp))

heatmap!(axQ, Qavg,          colormap = :balance)
heatmap!(axS, stats.Jˢ.avg,  colormap = :balance)
heatmap!(axu, stats.τx.avg,  colormap = :balance)
heatmap!(axv, stats.τy.avg,  colormap = :balance)

save("mean_fluxes.png", fig)

fig = Figure(resolution = (800, 600))
axQ = Axis(fig[1, 1], title = "Heat Flux")
axS = Axis(fig[1, 2], title = "Salinity Flux")
axu = Axis(fig[2, 1], title = "Zonal stress")
axv = Axis(fig[2, 2], title = "Meridional stress")

ρ  = reference_density(earth.model.ocean.model)
cp = heat_capacity(earth.model.ocean.model)

heatmap!(axQ, stats.Jᵀ.std * ρ * cp, colormap = :balance)
heatmap!(axS, stats.Jˢ.std,          colormap = :balance)
heatmap!(axu, stats.τx.std,          colormap = :balance)
heatmap!(axv, stats.τy.std,          colormap = :balance)

save("fluxes_std.png", fig)

fig = Figure(resolution = (800, 600))
axQ = Axis(fig[1, 1], title = "Heat Flux")
axS = Axis(fig[1, 2], title = "Salinity Flux")
axu = Axis(fig[2, 1], title = "Zonal stress")
axv = Axis(fig[2, 2], title = "Meridional stress")

ρ  = reference_density(earth.model.ocean.model)
cp = heat_capacity(earth.model.ocean.model)

heatmap!(axQ, stats.Jᵀ.max * ρ * cp, colormap = :balance)
heatmap!(axS, stats.Jˢ.max,          colormap = :balance)
heatmap!(axu, stats.τx.max,          colormap = :balance)
heatmap!(axv, stats.τy.max,          colormap = :balance)

save("fluxes_max.png", fig)

fig = Figure(resolution = (800, 600))
axQ = Axis(fig[1, 1], title = "Heat Flux")
axS = Axis(fig[1, 2], title = "Salinity Flux")
axu = Axis(fig[2, 1], title = "Zonal stress")
axv = Axis(fig[2, 2], title = "Meridional stress")

ρ  = reference_density(earth.model.ocean.model)
cp = heat_capacity(earth.model.ocean.model)

heatmap!(axQ, stats.Jᵀ.min * ρ * cp, colormap = :balance)
heatmap!(axS, stats.Jˢ.min,          colormap = :balance)
heatmap!(axu, stats.τx.min,          colormap = :balance)
heatmap!(axv, stats.τy.min,          colormap = :balance)

save("fluxes_min.png", fig)
