using ClimaOcean
using ClimaOcean.DataWrangling.JRA55: JRA55_field_time_series
using Oceananigans
using Tidejinks
using Dates
using GLMakie
import SPICE

backend = JRA55NetCDFBackend(41)
pa = JRA55_field_time_series(:sea_level_pressure; backend)

pa_with_tides = FieldTimeSeries{Center, Center, Nothing}(pa.grid, pa.times,
                                                         path = "sea_level_pressure_plus_tides.jld2",
                                                         name = "pa",
                                                         backend = OnDisk())

ρₒΦt = FieldTimeSeries{Center, Center, Nothing}(pa.grid, pa.times,
                                                path = "tidal_force_potential.jld2",
                                                name = "ρₒΦ",
                                                backend = OnDisk())

grid = pa.grid
ρₒ = 1020
Φ = Field{Center, Center, Nothing}(grid)
ρₒΦ = Field(ρₒ * Φ)

kernel_meta_file = "kernels.txt"
Tidejinks.wrangle_spice_kernels(kernel_meta_file)
SPICE.furnsh(kernel_meta_file)

times = pa.times
Nt = length(pa)

for n = 1:17 # Nt
    @show dt = Second(times[n])
    t = DateTime(1993, 1, 1, 1) + dt
    @info "Computing tides at $t"
    Tidejinks.compute_tidal_potential!(Φ, t)
    compute!(ρₒΦ)
    set!(ρₒΦt, ρₒΦ, n) #, pa.times[n])

    parent(ρₒΦ) .+= parent(pa[n])
    set!(pa_with_tides, ρₒΦ, n) #, pa.times[n])
end

using Statistics

fig = Figure()

ax1 = Axis(fig[1, 1], title="Pressure")
ax2 = Axis(fig[1, 2], title="Tides")

slider = Slider(fig[3, 1:2], startvalue=1, range=1:17)
n = slider.value #Observable(1)

titlestr = @lift string(pa.times[$n] ./ Oceananigans.Units.days)
Label(fig[0, 1:2], titlestr)

#pₐ_plus_ρₒΦ = @lift pa_with_tides[$n]

pₐ = @lift begin
    p = pa[$n]
    interior(p, :, :, 1) .- 101325 #mean(p)
end

mean_ρₒΦ = mean(ρₒΦt)

ρₒΦ = @lift begin
    ρₒΦ = ρₒΦt[$n]
    interior(ρₒΦ, :, :, 1) .- mean_ρₒΦ
end

#hm = heatmap!(ax1, pₐ_plus_ρₒΦ)
#Colorbar(fig[2, 1], hm, vertical=false)

hm = heatmap!(ax1, pₐ, colormap=:balance, colorrange=(-4000, 4000))
Colorbar(fig[2, 1], hm, vertical=false)

hm = heatmap!(ax2, ρₒΦ, colormap=:balance, colorrange=(-4000, 4000))
Colorbar(fig[2, 2], hm, vertical=false)

display(fig)

# heatmap(Φ)

