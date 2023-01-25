using Oceananigans
using Oceananigans.BuoyancyModels: ∂z_b, buoyancy_perturbation
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using JLD2
using GLMakie
using Printf

dir = "../data" #/storage1/greg"
# closure_name = "RiBasedVerticalDiffusivity"
closure_name = "CATKEVerticalDiffusivity"
#filename = "near_global_360_150_48_$(closure_name)_fields.jld2"
filename = "near_global_1440_600_48_$(closure_name)_fields.jld2"
filepath = joinpath(dir, filename)

file = jldopen(filepath)
reference_density = file["buoyancy/model/equation_of_state/reference_density"]
close(file)

equation_of_state = TEOS10EquationOfState(; reference_density)
buoyancy = SeawaterBuoyancy(; equation_of_state)

ut = FieldTimeSeries(filepath, "u")
vt = FieldTimeSeries(filepath, "v")
Tt = FieldTimeSeries(filepath, "T")
St = FieldTimeSeries(filepath, "S")

xu, yu, zu = nodes(ut)
xv, yv, zv = nodes(vt)
xc, yc, zc = nodes(Tt)
zw = znodes(Face, ut.grid)

times = Tt.times
grid = Tt.grid
Nx, Ny, Nz = size(grid)
Nt = length(times)

fig = Figure(resolution=(2400, 1800))

ax_T = Axis(fig[3, 2], xlabel="Longitude", ylabel="Latitude")
ax_S = Axis(fig[4, 2], xlabel="Longitude", ylabel="Latitude")
ax_u = Axis(fig[3, 3], xlabel="Longitude", ylabel="Latitude")
ax_v = Axis(fig[4, 3], xlabel="Longitude", ylabel="Latitude")

slider = Slider(fig[2, 2:3], range=1:Nt, startvalue=1)
n = slider.value

title = @lift string("Near-global one degree ocean simulation after ",
                     prettytime(times[$n] - times[1]))
Label(fig[1, :], title)

T = @lift interior(Tt[$n], :, :, 1)
S = @lift interior(St[$n], :, :, 1)
u = @lift interior(ut[$n], :, :, 1)
v = @lift interior(vt[$n], :, :, 1)

hm_T = heatmap!(ax_T, xc, yc, T, colormap=:thermal, colorrange=(2, 30))
hm_S = heatmap!(ax_S, xc, yc, S, colormap=:haline, colorrange=(30, 37))
hm_u = heatmap!(ax_u, xu, yu, u, colormap=:balance, colorrange=(-2e-1, 2e-1))
hm_v = heatmap!(ax_v, xv, yv, v, colormap=:balance, colorrange=(-2e-1, 2e-1))

Colorbar(fig[3, 1], hm_T; vertical=true, tellheight=false, flipaxis=false, label="Temperature (ᵒC)")
Colorbar(fig[4, 1], hm_S; vertical=true, tellheight=false, flipaxis=false, label="Salinity (g kg⁻¹)")
Colorbar(fig[3, 4], hm_u; vertical=true, tellheight=false, flipaxis=true, label="Zonal velocity (m s⁻¹)")
Colorbar(fig[4, 4], hm_v; vertical=true, tellheight=false, flipaxis=true, label="Meridional velocity (m s⁻¹)")

display(fig)

# record(fig, "one_degree_near_global_simulation_$(closure_name).mp4", 1:Nt, framerate=24) do nn
#     n[] = nn
# end

