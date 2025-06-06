px_per_unit = 4

ρ  = reference_density(earth.model.ocean)
cp = heat_capacity(earth.model.ocean)

using GLMakie

fig = Figure(size = (1200, 900))

axQ  = Axis(fig[1, 1], title = "Heat Flux")
axS  = Axis(fig[1, 2], title = "Salinity Flux")
axQc = Axis(fig[2, 1], title = "Sensible Heat Flux")
axQv = Axis(fig[2, 2], title = "Latent Heat Flux")
axu  = Axis(fig[3, 1], title = "Zonal stress")
axv  = Axis(fig[3, 2], title = "Meridional stress")

hmQ  = heatmap!(axQ,  stats.Jᵀ.mean * ρ * cp, colormap = :balance, colorrange = (-100, 250))
hmS  = heatmap!(axS,  stats.Jˢ.mean,          colormap = :balance, colorrange = (-2.5e-6, 2.5e-6))
hmQc = heatmap!(axQc, stats.Qc.mean,          colormap = :balance, colorrange = (-100, 250))
hmQv = heatmap!(axQv, stats.Qv.mean,          colormap = :balance, colorrange = (-100, 250))
hmu  = heatmap!(axu,  stats.τx.mean * ρ,      colormap = :balance, colorrange = (-0.2, 0.2))
hmv  = heatmap!(axv,  stats.τy.mean * ρ,      colormap = :balance, colorrange = (-0.2, 0.2))

Colorbar(fig[1, 0], hmQ,  flipaxis=false, label = "W/m²")
Colorbar(fig[1, 3], hmS,  flipaxis=true,  label = "kg/m²/s")
Colorbar(fig[2, 0], hmQc, flipaxis=false, label = "W/m²")
Colorbar(fig[2, 3], hmQv, flipaxis=true,  label = "W/m²")
Colorbar(fig[3, 0], hmu,  flipaxis=false, label = "N/m²")
Colorbar(fig[3, 3], hmv,  flipaxis=true,  label = "N/m²")

save("fluxes_mean.png", fig; px_per_unit)


fig = Figure(size = (1200, 900))

axQ  = Axis(fig[1, 1], title = "Heat Flux")
axS  = Axis(fig[1, 2], title = "Salinity Flux")
axQc = Axis(fig[2, 1], title = "Sensible Heat Flux")
axQv = Axis(fig[2, 2], title = "Latent Heat Flux")
axu  = Axis(fig[3, 1], title = "Zonal stress")
axv  = Axis(fig[3, 2], title = "Meridional stress")

hmQ  = heatmap!(axQ,  stats.Jᵀ.std * ρ * cp, colormap = :deep, colorrange = (0, 200))
hmS  = heatmap!(axS,  stats.Jˢ.std,          colormap = :deep, colorrange = (0, 5e-6))
hmQc = heatmap!(axQc, stats.Qc.std,          colormap = :deep, colorrange = (0, 200))
hmQv = heatmap!(axQv, stats.Qv.std,          colormap = :deep, colorrange = (0, 200))
hmu  = heatmap!(axu,  stats.τx.std * ρ,      colormap = :deep, colorrange = (0, 0.3))
hmv  = heatmap!(axv,  stats.τy.std * ρ,      colormap = :deep, colorrange = (0, 0.3))

Colorbar(fig[1, 0], hmQ,  flipaxis=false, label = "W/m²")
Colorbar(fig[1, 3], hmS,  flipaxis=true,  label = "kg/m²/s")
Colorbar(fig[2, 0], hmQc, flipaxis=false, label = "W/m²")
Colorbar(fig[2, 3], hmQv, flipaxis=true,  label = "W/m²")
Colorbar(fig[3, 0], hmu,  flipaxis=false, label = "N/m²")
Colorbar(fig[3, 3], hmv,  flipaxis=true,  label = "N/m²")

save("fluxes_std.png", fig; px_per_unit)


fig = Figure(size = (1200, 900))

axQ  = Axis(fig[1, 1], title = "Heat Flux")
axS  = Axis(fig[1, 2], title = "Salinity Flux")
axQc = Axis(fig[2, 1], title = "Sensible Heat Flux")
axQv = Axis(fig[2, 2], title = "Latent Heat Flux")
axu  = Axis(fig[3, 1], title = "Zonal stress")
axv  = Axis(fig[3, 2], title = "Meridional stress")

hmQ  = heatmap!(axQ,  stats.Jᵀ.max * ρ * cp, colormap = :magma, colorrange = (0,  1000))
hmS  = heatmap!(axS,  stats.Jˢ.max,          colormap = :magma, colorrange = (1e-6, 5e-5))
hmQc = heatmap!(axQc, stats.Qc.max,          colormap = :magma, colorrange = (0, 1000))
hmQv = heatmap!(axQv, stats.Qv.max,          colormap = :magma, colorrange = (0, 1000))
hmu  = heatmap!(axu,  stats.τx.max * ρ,      colormap = :magma, colorrange = (0.1,  2.0))
hmv  = heatmap!(axv,  stats.τy.max * ρ,      colormap = :magma, colorrange = (0.1,  2.0))

Colorbar(fig[1, 0], hmQ,  flipaxis=false, label = "W/m²")
Colorbar(fig[1, 3], hmS,  flipaxis=true,  label = "kg/m²/s")
Colorbar(fig[2, 0], hmQc, flipaxis=false, label = "W/m²")
Colorbar(fig[2, 3], hmQv, flipaxis=true,  label = "W/m²")
Colorbar(fig[3, 0], hmu,  flipaxis=false, label = "N/m²")
Colorbar(fig[3, 3], hmv,  flipaxis=true,  label = "N/m²")

save("fluxes_max.png", fig; px_per_unit)


fig = Figure(size = (1200, 900))

axQ  = Axis(fig[1, 1], title = "Heat Flux")
axS  = Axis(fig[1, 2], title = "Salinity Flux")
axQc = Axis(fig[2, 1], title = "Sensible Heat Flux")
axQv = Axis(fig[2, 2], title = "Latent Heat Flux")
axu  = Axis(fig[3, 1], title = "Zonal stress")
axv  = Axis(fig[3, 2], title = "Meridional stress")

hmQ  = heatmap!(axQ,  stats.Jᵀ.min * ρ * cp, colormap = :haline, colorrange = (-1500,  0))
hmS  = heatmap!(axS,  stats.Jˢ.min,          colormap = :haline, colorrange = (-5e-6,  0))
hmQc = heatmap!(axQc, stats.Qc.min,          colormap = :haline, colorrange = (-1000, 0))
hmQv = heatmap!(axQv, stats.Qv.min,          colormap = :haline, colorrange = (-1000, 0))
hmu  = heatmap!(axu,  stats.τx.min * ρ,      colormap = :haline, colorrange = (-2.0,  0.1))
hmv  = heatmap!(axv,  stats.τy.min * ρ,      colormap = :haline, colorrange = (-2.0,  0.1))

Colorbar(fig[1, 0], hmQ,  flipaxis=false, label = "W/m²")
Colorbar(fig[1, 3], hmS,  flipaxis=true,  label = "kg/m²/s")
Colorbar(fig[2, 0], hmQc, flipaxis=false, label = "W/m²")
Colorbar(fig[2, 3], hmQv, flipaxis=false,  label = "W/m²")
Colorbar(fig[3, 0], hmu,  flipaxis=false, label = "N/m²")
Colorbar(fig[3, 3], hmv,  flipaxis=true,  label = "N/m²")

save("fluxes_min.png", fig; px_per_unit)
