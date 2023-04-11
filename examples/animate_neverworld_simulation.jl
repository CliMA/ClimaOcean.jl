using Oceananigans
using Oceananigans.Units
using GLMakie
using Printf
using Statistics

# Some plotting parameters
ζlim = 6e-5

# Plot a limited number of years (here, 10)
Nyears = 3
t_end = 200 * 360 * day # 200 "years"
t_start = (200 - Nyears) * 360 * day # 200 "years"

# Load OnDisk time-series without loading massive amount of data
backend = OnDisk()
bxyt = FieldTimeSeries("neverworld_xy.jld2", "b"; backend)

# Determine which iterations to load into memory
t̃ = bxyt.times
n_end   = searchsortedfirst(t̃, t_end)
n_start = searchsortedfirst(t̃, t_start) - 1

t = times = t̃[n_start:n_end]
Nt = length(t)

bxyt = FieldTimeSeries("neverworld_xy.jld2", "b"; times)
uxyt = FieldTimeSeries("neverworld_xy.jld2", "u"; times)
vxyt = FieldTimeSeries("neverworld_xy.jld2", "v"; times)
exyt = FieldTimeSeries("neverworld_xy.jld2", "e"; times)
κxyt = FieldTimeSeries("neverworld_xy.jld2", "κᶜ"; times)

byzt = FieldTimeSeries("neverworld_zonal_average.jld2", "b"; times)
eyzt = FieldTimeSeries("neverworld_zonal_average.jld2", "e"; times)
κyzt = FieldTimeSeries("neverworld_zonal_average.jld2", "κᶜ"; times)

# Hack to set immersed regions to NaN
for ft in (bxyt, uxyt, exyt, κxyt, byzt, eyzt, κyzt)
    fp = parent(ft)
    fp[fp .== 0] .= NaN
end

fig = Figure(resolution=(2000, 1800))

axbxy = Axis(fig[2, 1], xlabel="Longitude", ylabel="Latitude", title="b(x, y)")
axζxy = Axis(fig[2, 2], xlabel="Longitude", ylabel="Latitude", title="ζ(x, y)")
axexy = Axis(fig[2, 3], xlabel="Longitude", ylabel="Latitude", title="e(x, y)")

axNyz = Axis(fig[3, 1], xlabel="Longitude", ylabel="z (m)", title="N²(y, z)")
axeyz = Axis(fig[3, 2], xlabel="Longitude", ylabel="z (m)", title="e(y, z)")
axκyz = Axis(fig[3, 3], xlabel="Longitude", ylabel="z (m)", title="κᶜ(y, z)")

slider = Slider(fig[5, 1:3], range=1:Nt, startvalue=Nt)
n = slider.value

title = @lift begin
    t_day = t[$n] / day
    t_year = floor(Int, t_day / 360)
    yearday = t_day % 360
    return @sprintf("CATKE Neverworld at % 3d years, % 3d days", t_year, yearday)
end

Label(fig[0, 1:3], title, fontsize=36)

using Oceananigans.Operators: ζ₃ᶠᶠᶜ
grid = bxyt.grid
grid = grid.underlying_grid
Nx, Ny, Nz = size(grid)
uxy = Field{Face, Center, Center}(grid, indices=(:, :, 1))
vxy = Field{Center, Face, Center}(grid, indices=(:, :, 1))
ζop = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, uxy, vxy)
ζxy = Field(ζop, indices=(:, :, 1))

ζxyn = @lift begin
    uxyn = uxyt[$n]
    vxyn = vxyt[$n]
    interior(uxy) .= interior(uxyn)
    interior(vxy) .= interior(vxyn)
    compute!(ζxy)
    return interior(ζxy, :, :, 1)
end

byz = Field{Nothing, Center, Center}(grid)
N²yz = Field(∂z(byz))

N²yzn = @lift begin
    byzn = byzt[$n]
    interior(byz) .= interior(byzn)
    compute!(N²yz)
    return interior(N²yz, 1, :, :)
end

bxyn = @lift interior(bxyt[$n], :, :, 1)
exyn = @lift interior(exyt[$n], :, :, 1)

byzn = @lift interior(byzt[$n], 1, :, :)
eyzn = @lift interior(eyzt[$n], 1, :, :)
κyzn = @lift interior(κyzt[$n], 1, :, :)

κlims = @lift begin
    κ_mean = mean(filter(!isnan, interior(κyzt[$n])))
    (κ_mean/2, 3κ_mean)
end

eyzlims = @lift begin
    e_mean = mean(filter(!isnan, interior(eyzt[$n])))
    (e_mean/2, 3e_mean)
end

exylims = @lift begin
    e_mean = mean(filter(!isnan, interior(eyzt[$n])))
    (e_mean/10, e_mean)
end

x, y, z = nodes(bxyt)
xf, yf, zf = nodes(grid, Face(), Face(), Face())
cbkw = (vertical=false, flipaxis=true, tellwidth=false)

hm = heatmap!(axbxy, x, y, bxyn, colorrange=(-0.05, 0.05))
Colorbar(fig[1, 1], hm; label="Buoyancy", cbkw...)

hm = heatmap!(axζxy, x, y, ζxyn, colormap=:balance, colorrange=(-ζlim, ζlim))
Colorbar(fig[1, 2], hm; label="Relative vertical vorticity (s⁻¹)", cbkw...)

hm = heatmap!(axexy, x, y, exyn, colormap=:solar, colorrange=exylims)
Colorbar(fig[1, 3], hm; label="Turbulent kinetic energy (m² s⁻²)", cbkw...)

cbkw = (vertical=false, flipaxis=false, tellwidth=false)

hm = heatmap!(axNyz, y, zf, N²yzn, nan_color=:gray, colorrange=(1e-6, 1e-4))
Colorbar(fig[4, 1], hm; label="Buoyancy gradient (s⁻²)", cbkw...)
contour!(axNyz, y, z, byzn, levels=30, color=:white, linewidth=3)

hm = heatmap!(axeyz, y, z, eyzn, nan_color=:gray, colormap=:solar, colorrange=eyzlims)
Colorbar(fig[4, 2], hm; label="Turbulent kinetic energy (m² s⁻²)", cbkw...)
contour!(axeyz, y, z, byzn, levels=30, color=:white, linewidth=3)

hm = heatmap!(axκyz, y, zf, κyzn, nan_color=:gray, colormap=:solar, colorrange=κlims)
Colorbar(fig[4, 3], hm; label="Tracer eddy diffusivity (m² s⁻¹)", cbkw...)
contour!(axκyz, y, z, byzn, levels=30, color=:white, linewidth=3)

ylims!(axNyz, -4000, 0)
ylims!(axeyz, -4000, 0)
ylims!(axκyz, -4000, 0)

display(fig)

record(fig, "catke_neverworld.mp4", 1:Nt, framerate=8) do nn
     @info "Recording frame $nn of $Nt..."
     n[] = nn
end

