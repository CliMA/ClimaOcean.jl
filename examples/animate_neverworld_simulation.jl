using Oceananigans
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.Units
using GLMakie
using Printf

#=
bxyt = FieldTimeSeries("neverworld_xy.jld2", "b")
uxyt = FieldTimeSeries("neverworld_xy.jld2", "u")
ζxyt = FieldTimeSeries("neverworld_xy.jld2", "ζ")
exyt = FieldTimeSeries("neverworld_xy.jld2", "e")
κxyt = FieldTimeSeries("neverworld_xy.jld2", "κᶜ")

byzt = FieldTimeSeries("neverworld_yz.jld2", "b")
uyzt = FieldTimeSeries("neverworld_yz.jld2", "u")
ζyzt = FieldTimeSeries("neverworld_yz.jld2", "ζ")
eyzt = FieldTimeSeries("neverworld_yz.jld2", "e")
κyzt = FieldTimeSeries("neverworld_yz.jld2", "κᶜ")

t = bxyt.times
Nt = length(t)

for ft in (bxyt, uxyt, ζxyt, exyt, κxyt, byzt, uyzt, ζyzt, eyzt, κyzt)
    fp = parent(ft)
    fp[fp .== 0] .= NaN
end
=#

fig = Figure(resolution=(3200, 1800))

axbxy = Axis(fig[2, 1], xlabel="Longitude", ylabel="Latitude", title="b(x, y)")
axuxy = Axis(fig[2, 2], xlabel="Longitude", ylabel="Latitude", title="u(x, y)")
axζxy = Axis(fig[2, 3], xlabel="Longitude", ylabel="Latitude", title="ζ(x, y)")
axexy = Axis(fig[2, 4], xlabel="Longitude", ylabel="Latitude", title="e(x, y)")
axκxy = Axis(fig[2, 5], xlabel="Longitude", ylabel="Latitude", title="κᶜ(x, y)")

axbyz = Axis(fig[3, 1], xlabel="Longitude", ylabel="z (m)", title="b(y, z)")
axuyz = Axis(fig[3, 2], xlabel="Longitude", ylabel="z (m)", title="u(y, z)")
axζyz = Axis(fig[3, 3], xlabel="Longitude", ylabel="z (m)", title="ζ(y, z)")
axeyz = Axis(fig[3, 4], xlabel="Longitude", ylabel="z (m)", title="e(y, z)")
axκyz = Axis(fig[3, 5], xlabel="Longitude", ylabel="z (m)", title="κᶜ(y, z)")

slider = Slider(fig[4, 1:5], range=1:Nt, startvalue=Nt)
n = slider.value

title = @lift @sprintf("CATKE Neverworld at t = % 3d days", t[$n] / day)
Label(fig[0, 1:5], title, fontsize=36)

clip(x) = max(zero(x), x)
fᵉ = identity #log10 ∘ clip
fᵏ = identity #log10 ∘ clip

bxyn = @lift interior(bxyt[$n], :, :, 1)
uxyn = @lift interior(uxyt[$n], :, :, 1)
ζxyn = @lift interior(ζxyt[$n], :, :, 1)
exyn = @lift fᵉ.(interior(exyt[$n], :, :, 1))
κxyn = @lift fᵏ.(interior(κxyt[$n], :, :, 1))

byzn = @lift interior(byzt[$n], 1, :, :)
uyzn = @lift interior(uyzt[$n], 1, :, :)
ζyzn = @lift interior(ζyzt[$n], 1, :, :)
eyzn = @lift fᵉ.(interior(eyzt[$n], 1, :, :))
κyzn = @lift fᵏ.(interior(κyzt[$n], 1, :, :))

ulim = 2.0
ζlim = 6e-5
κlims = (1e-2, 2e2)
elims = (1e-4, 6e-3)
x, y, z = nodes(bxyt)
cbkw = (vertical=false, flipaxis=true, tellwidth=false)

hm = heatmap!(axbxy, x, y, bxyn, colorrange=(-0.02, 0.06))
Colorbar(fig[1, 1], hm; label="Buoyancy", cbkw...)

hm = heatmap!(axuxy, x, y, uxyn, colormap=:balance, colorrange=(-ulim, ulim))
Colorbar(fig[1, 2], hm; label="Zonal velocity (m s⁻¹)", cbkw...)

hm = heatmap!(axζxy, x, y, ζxyn, colormap=:balance, colorrange=(-ζlim, ζlim))
Colorbar(fig[1, 3], hm; label="Relative vertical vorticity (s⁻¹)", cbkw...)

hm = heatmap!(axexy, x, y, exyn, colormap=:solar, colorrange=elims)
Colorbar(fig[1, 4], hm; label="Turbulent kinetic energy (m² s⁻²)", cbkw...)

hm = heatmap!(axκxy, x, y, κxyn, colormap=:solar, colorrange=κlims)
Colorbar(fig[1, 5], hm; label="Vertical tracer diffusivity (m² s⁻¹)", cbkw...)

hm = heatmap!(axbyz, y, z, byzn, nan_color=:gray, colorrange=(-0.02, 0.06))
hm = heatmap!(axuyz, y, z, uyzn, nan_color=:gray, colormap=:balance, colorrange=(-ulim, ulim))
hm = heatmap!(axζyz, y, z, ζyzn, nan_color=:gray, colormap=:balance, colorrange=(-ζlim, ζlim))
hm = heatmap!(axeyz, y, z, eyzn, nan_color=:gray, colormap=:solar, colorrange=elims)
hm = heatmap!(axκyz, y, z, κyzn, nan_color=:gray, colormap=:solar, colorrange=κlims)

display(fig)

record(fig, "catke_neverworld.mp4", 1:Nt, framerate=36) do nn
    n[] = nn
end

