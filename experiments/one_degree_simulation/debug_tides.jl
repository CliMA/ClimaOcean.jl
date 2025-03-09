using ClimaOcean
using ClimaOcean.DataWrangling.JRA55: JRA55_field_time_series
using Oceananigans
using Tidejinks
using Dates
using GLMakie
import SPICE
using Statistics

backend = JRA55NetCDFBackend(41)
ρₒ = 1020
pa = JRA55_field_time_series(:sea_level_pressure; backend)
grid = pa.grid

# dt = 30 * 60 # minutes
# times = 0.0:dt:(30 * 86400.0)
times = pa.times

Φ☾t = FieldTimeSeries{Center, Center, Nothing}(pa.grid, times)
Φpt = FieldTimeSeries{Center, Center, Nothing}(pa.grid, times)

grid = pa.grid
Φ = Field{Center, Center, Nothing}(grid)

kernel_meta_file = "kernels.txt"
Tidejinks.wrangle_spice_kernels(kernel_meta_file)
SPICE.furnsh(kernel_meta_file)

Nt = length(times)

for n = 1:Nt
    Δt = Second(times[n])
    t = DateTime(1993, 1, 1, 1) + Δt
    @info "Computing tides at $t"
    Tidejinks.compute_tidal_potential!(Φ, t)
    parent(Φ☾t[n]) .= parent(Φ)

    pan = parent(pa[n])
    parent(Φpt[n]) .= (pan .- mean(pan)) ./ ρₒ
end

#=
fig = Figure()
ax = Axis(fig[1, 1])
slider = Slider(fig[2, 1], startvalue=1, range=1:Nt)
n = slider.value #Observable(1)

titlestr = @lift string(times[$n] ./ Oceananigans.Units.days)
Label(fig[0, 1], titlestr, tellwidth=false)
Φn = @lift Φt[$n]
hm = heatmap!(ax, Φn)
display(fig)
=#

# heatmap(Φ)

