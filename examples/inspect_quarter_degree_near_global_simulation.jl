using Oceananigans
using Oceananigans.BuoyancyModels: ∂z_b, buoyancy_perturbation
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using JLD2
using GLMakie
using Printf

#filepath = "near_global_1440_600_87_RiBasedVerticalDiffusivity_fields_surface.jld2"
#filepath = "catke_test_near_global_1440_600_87_fields_surface.jld2"
#filepath = "catke_flat_bottom_fields_surface.jld2"
filepath = "ri_based_fields_surface.jld2"

file = jldopen(filepath)
@show keys(file["buoyancy/model/equation_of_state"])
close(file)

#=
equation_of_state = TEOS10EquationOfState(; reference_density)
buoyancy = SeawaterBuoyancy(; equation_of_state)
=#

ut = FieldTimeSeries(filepath, "u")
vt = FieldTimeSeries(filepath, "v")
Tt = FieldTimeSeries(filepath, "T")
St = FieldTimeSeries(filepath, "S")
ζt = FieldTimeSeries(filepath, "ζ")
et = FieldTimeSeries(filepath, "T")

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
ax_e = Axis(fig[5, 2], xlabel="Longitude", ylabel="Latitude")
ax_z = Axis(fig[5, 3], xlabel="Longitude", ylabel="Latitude")

slider = Slider(fig[2, 2:3], range=1:Nt, startvalue=1)
n = slider.value

title = @lift string("Near-global quarter degree ocean simulation after ",
                     prettytime(times[$n] - times[1]))
Label(fig[1, :], title)

Nx, Ny, Nz = size(grid)
land = findall(h -> h > 0, grid.immersed_boundary.bottom_height[1:Nx, 1:Ny])

for n = 1:Nt
    en = interior(et[n], :, :, 1)
    en[land] .= NaN
end

T = @lift interior(Tt[$n], :, :, 1)
S = @lift interior(St[$n], :, :, 1)
u = @lift interior(ut[$n], :, :, 1)
v = @lift interior(vt[$n], :, :, 1)
e = @lift interior(et[$n], :, :, 1)
ζ = @lift interior(ζt[$n], :, :, 1)

hm_T = heatmap!(ax_T, xc, yc, T, colormap=:thermal, colorrange=(2, 30))
hm_S = heatmap!(ax_S, xc, yc, S, colormap=:haline,  colorrange=(2, 32))
hm_u = heatmap!(ax_u, xu, yu, u, colormap=:balance, colorrange=(-2e-1, 2e-1))
hm_v = heatmap!(ax_v, xv, yv, v, colormap=:balance, colorrange=(-2e-1, 2e-1))
hm_e = heatmap!(ax_e, xc, yc, e, colormap=:solar,   colorrange=(0, 10.0))
hm_z = heatmap!(ax_z, xu, yv, ζ, colormap=:balance, colorrange=(-5e-5, 5e-5))

Colorbar(fig[3, 1], hm_T; vertical=true, tellheight=false, flipaxis=false, label="Temperature (ᵒC)")
Colorbar(fig[4, 1], hm_S; vertical=true, tellheight=false, flipaxis=false, label="Salinity (g kg⁻¹)")
Colorbar(fig[3, 4], hm_u; vertical=true, tellheight=false, flipaxis=true,  label="Zonal velocity (m s⁻¹)")
Colorbar(fig[4, 4], hm_v; vertical=true, tellheight=false, flipaxis=true,  label="Meridional velocity (m s⁻¹)")
Colorbar(fig[5, 1], hm_e; vertical=true, tellheight=false, flipaxis=false, label="TKE (m² s⁻²)")
Colorbar(fig[5, 4], hm_z; vertical=true, tellheight=false, flipaxis=true,  label="Vorticity (s⁻¹)")

display(fig)

# record(fig, "one_degree_near_global_simulation_$(closure_name).mp4", 1:Nt, framerate=24) do nn
#     n[] = nn
# end

