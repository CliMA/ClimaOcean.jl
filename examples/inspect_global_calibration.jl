using Oceananigans
using Oceananigans.BuoyancyModels: ∂z_b, buoyancy_perturbation
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using JLD2
using GLMakie
using Printf
using Statistics
using ColorSchemes

set_theme!(Theme(fontsize=36))

Nens = 6
Niters = 10
dir = "../data" #/storage1/greg"
prefix = "gm_one_degree_calibration_depth_1007_meters_eki_iteration"

function filepath(iter, particle)
    if iter == 0
        num = 0
    else
        num = iter * Nens + particle - 6
    end

    filename = string(prefix, num, "_particle", particle, ".jld2")

    return joinpath(dir, filename)
end

file = jldopen(filepath(1, 1))
@show κ_skew = file["closure"]["4"]["κ_skew"]
@show κ_symmetric = file["closure"]["4"]["κ_symmetric"]
close(file)

Tt = FieldTimeSeries(filepath(1, 1), "T")
Tref = Tt[2]
grid = Tt.grid
x, y, z = nodes(Tt)

Nx, Ny, Nz = size(grid)

Tens = zeros(Nx, Ny, Nens)
Tstd = zeros(Nx, Ny, Niters)

κG = [zeros(Nens) for i = 1:Niters]
κR = [zeros(Nens) for i = 1:Niters]

for iter = 1:Niters
    Tens .= 0

    for k = 1:6
        T = FieldTimeSeries(filepath(iter, k), "T")[2]
        Tens[:, :, k] .= interior(T, :, :, 1)

        file = jldopen(filepath(iter, k))
        κG[iter][k] = file["closure"]["4"]["κ_skew"]
        κR[iter][k] = file["closure"]["4"]["κ_symmetric"]
        close(file)
    end

    Tstd[:, :, iter] .= std(Tens, dims=3)[:, :, 1]

    Ti = view(Tstd, :, :, iter:iter)
    Ti[interior(Tref) .== 0] .= NaN
end

Nplot = 5
ticks = [2000, 3000, 4000]

fig = Figure(resolution=(1800, 900))
axT = Axis(fig[2, 1:3], xlabel="Longitude", ylabel="Latitude")
axθ = Axis(fig[2, 4], xlabel="κᴳᴹ (m² s⁻¹)", ylabel="κᴿᵉᵈⁱ (m² s⁻¹)", aspect=1,
           xticks=ticks, yticks=ticks, yaxisposition=:right)

xlims!(axθ, 1000, 4500)
ylims!(axθ, 1000, 4500)

hidedecorations!(axT)

#colors = ColorSchemes.Paired_11.colors
colors = ColorSchemes.gist_earth.colors

for iter = 1:Nplot
    scatter!(axθ, κG[iter], κR[iter], color=(colors[iter], 0.4), markersize=15)
end

#slider = Slider(fig[3, 1:5], range=1:Nplot, startvalue=1)
#n = slider.value
n = Observable(1)

color_n = @lift colors[($n - 1) * 20 + 1]
κG_n = @lift κG[$n]
κR_n = @lift κR[$n]
scatter!(axθ, κG_n, κR_n, color=(color_n, 1.0), markersize=20)

Tlim = maximum(abs, filter(!isnan, Tstd)) / 10
T = @lift Tstd[:, :, $n]
hm = heatmap!(axT, x, y, T, colormap=:thermal, colorrange=(0, Tlim), nan_color=:black)

label = @lift string("Standard deviation of Tₖ(z=-1000 m) (ᵒC) at EKI iteration ", $n)
Colorbar(fig[1, 1:3], hm; vertical=false, flipaxis=true, label, tellheight=true,
         width=Relative(0.9))

display(fig)

record(fig, "global_calibration.gif", 1:Nplot, framerate=2) do nn
    n[] = nn
end

