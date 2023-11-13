using Oceananigans
using Oceananigans.BuoyancyModels: ∂z_b, buoyancy_perturbation
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using JLD2
using GLMakie
using Printf

dir = "../data" #/storage1/greg"
# closure_name = "RiBasedVerticalDiffusivity"
closure_name = "CATKEVerticalDiffusivity"
filename = "near_global_360_150_48_$(closure_name)_fields.jld2"
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

ax_T = Axis(fig[2, 2], xlabel="Longitude", ylabel="Latitude")
ax_S = Axis(fig[3, 2], xlabel="Longitude", ylabel="Latitude")
ax_u = Axis(fig[2, 3], xlabel="Longitude", ylabel="Latitude")
ax_v = Axis(fig[3, 3], xlabel="Longitude", ylabel="Latitude")

ax_N1 = Axis(fig[4, 2], xlabel="Latitude", ylabel="z (m)")
ax_N2 = Axis(fig[4, 3], xlabel="Latitude", ylabel="z (m)")

slider = Slider(fig[5, 2:3], range=1:Nt, startvalue=1)
n = slider.value

title = @lift string("Near-global one degree ocean simulation after ",
                     prettytime(times[$n] - times[1]))
Label(fig[1, 1:2], title)

T = @lift interior(Tt[$n], :, :, grid.Nz)
S = @lift interior(St[$n], :, :, grid.Nz)
u = @lift interior(ut[$n], :, :, grid.Nz)
v = @lift interior(vt[$n], :, :, grid.Nz)

function N²n(n, indices)
    T = Tt[n]
    S = St[n]
    tracers = (; T, S)
    N²_op = KernelFunctionOperation{Center, Center, Face}(∂z_b, grid,
                                                          computed_dependencies=(buoyancy, tracers))
    interior(compute!(Field(N²_op; indices)), 1, :, :)
end


function bn(n, indices)
    T = Tt[n]
    S = St[n]
    tracers = (; T, S)
    b_op = KernelFunctionOperation{Center, Center, Center}(buoyancy_perturbation, grid,
                                                           computed_dependencies=(buoyancy, tracers))
    interior(compute!(Field(b_op; indices)), 1, :, :)
end

k0 = 12
i1 = 11
i2 = 153
T1 = @lift interior(Tt[$n], i1, :, k0:Nz)
T2 = @lift interior(Tt[$n], i2, :, k0:Nz)

N²1 = @lift N²n($n, (i1, :, k0:Nz))
N²2 = @lift N²n($n, (i2, :, k0:Nz))

b1 = @lift bn($n, (i1, :, k0:Nz))
b2 = @lift bn($n, (i2, :, k0:Nz))

hm_T = heatmap!(ax_T, xc, yc, T, colormap=:thermal, colorrange=(2, 30))
hm_S = heatmap!(ax_S, xc, yc, S, colormap=:haline, colorrange=(30, 37))
hm_u = heatmap!(ax_u, xu, yu, u, colormap=:balance, colorrange=(-2e-1, 2e-1))
hm_v = heatmap!(ax_v, xv, yv, v, colormap=:balance, colorrange=(-2e-1, 2e-1))

vlines!(ax_T, xc[i1], color=:white)
vlines!(ax_T, xc[i2], color=:white)
vlines!(ax_S, xc[i1], color=:white)
vlines!(ax_S, xc[i2], color=:white)

hm_N1 = heatmap!(ax_N1, yc, zw[k0:Nz], N²1; colormap=:thermal, colorrange=(1e-6, 2e-5))
hm_N2 = heatmap!(ax_N2, yc, zw[k0:Nz], N²2; colormap=:thermal, colorrange=(1e-6, 2e-5))

b1₀ = bn(1, (i1, :, k0:Nz))
bmax = maximum(b1₀)
bmin = minimum(b1₀)
levels = range(bmin, stop=bmax, length=40)
contour!(ax_N1, yc, zc[k0:Nz], b1; levels, color=:white)

b2₀ = bn(1, (i2, :, k0:Nz))
bmax = maximum(b2₀)
bmin = minimum(b2₀)
levels = range(bmin, stop=bmax, length=40)
contour!(ax_N2, yc, zc[k0:Nz], b2; levels, color=:white)

λ1 = xc[i1]
λ2 = xc[i2]
N1_label = @sprintf("Buoyancy frequency (s⁻²) at λ = %.1fᵒ", λ1)
N2_label = @sprintf("Buoyancy frequency (s⁻²) at λ = %.1fᵒ", λ2)

Colorbar(fig[2, 1], hm_T; vertical=true, tellheight=false, flipaxis=false, label="Temperature (ᵒC)")
Colorbar(fig[3, 1], hm_S; vertical=true, tellheight=false, flipaxis=false, label="Salinity (g kg⁻¹)")
Colorbar(fig[4, 1], hm_N1; vertical=true, tellheight=false, flipaxis=false, label=N1_label)

Colorbar(fig[2, 4], hm_u; vertical=true, tellheight=false, flipaxis=true, label="Zonal velocity (m s⁻¹)")
Colorbar(fig[3, 4], hm_v; vertical=true, tellheight=false, flipaxis=true, label="Meridional velocity (m s⁻¹)")
Colorbar(fig[4, 4], hm_N2; vertical=true, tellheight=false, flipaxis=true, label=N2_label)

display(fig)

record(fig, "one_degree_near_global_simulation_$(closure_name).mp4", 1:Nt, framerate=24) do nn
    n[] = nn
end

