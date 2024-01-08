using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels: buoyancy_frequency

using ClimaOcean
using ClimaOcean.DataWrangling.JRA55: jra55_prescribed_atmosphere
using ClimaOcean.DataWrangling.ECCO2: ecco2_field

using GLMakie
using Printf
using Dates

start_time = time_ns()

include("single_column_omip_ocean_component.jl")

epoch = Date(1992, 1, 1)
date = Date(1992, 10, 01)
start_seconds = Second(date - epoch).value
uᵢ = ecco2_field(:u_velocity, date)
vᵢ = ecco2_field(:v_velocity, date)
Tᵢ = ecco2_field(:temperature, date)
Sᵢ = ecco2_field(:salinity, date)

land = interior(Tᵢ) .< -10
interior(Tᵢ)[land] .= NaN
interior(Sᵢ)[land] .= NaN

teos10 = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=teos10)
tracers = (T=Tᵢ, S=Sᵢ)
N²_op = buoyancy_frequency(buoyancy, Tᵢ.grid, tracers)
N² = Field(N²_op)
compute!(N²)

elapsed = time_ns() - start_time
@info "Initial condition built. " * prettytime(elapsed * 1e-9)
start_time = time_ns()

#####
##### Construct the grid
#####

zc = znodes(Tᵢ)
zf = znodes(N²)

arch = CPU()

Δ = 1/4 # resolution in degrees
φ₁ = -90 + Δ/2
φ₂ = +90 - Δ/2
λ₁ = 0   + Δ/2
λ₂ = 360 - Δ/2
φe = φ₁:Δ:φ₂
λe = λ₁:Δ:λ₂

Nz = size(Tᵢ, 3)
fig = Figure(resolution=(1200, 1200))
map = Axis(fig[1, 1:3], xlabel="λ (degrees)", ylabel="φ (degrees)")
hm = heatmap!(map, λe, φe, interior(Tᵢ, :, :, Nz), colorrange=(0, 30), nan_color=:gray)
Colorbar(fig[1, 4], hm, label="Surface temperature (ᵒC)")

axT = Axis(fig[2, 1], ylabel="z (m)", xlabel="Temperature (ᵒC)")
axS = Axis(fig[2, 2], ylabel="z (m)", xlabel="Salinity (g/kg)")
axN = Axis(fig[2, 3], ylabel="z (m)", xlabel="Buoyancy frequency (s⁻²)")

φs = [50,   55, 0,   -30, -65, 34]
λs = [215, 310, 210, 160, 160, 34]

# Mediterranean locations
# λs = [34, 33, 5, 20, 30]
# φs = [34, 32, 38, 35, 33]

Nc = length(φs)

for n = 1:Nc
    local φ★
    local λ★
    local i★
    local j★

    φ★ = φs[n]
    λ★ = λs[n]

    i★ = searchsortedfirst(λe, λ★)
    j★ = searchsortedfirst(φe, φ★)

    scatter!(map, λ★, φ★, strokewidth=4, strokecolor=:black,
             color=:pink, markersize=20)

    label = string("λ = ", λ★, ", φ = ", φ★)
    scatterlines!(axT, interior(Tᵢ, i★, j★, :), zc; label)
    scatterlines!(axS, interior(Sᵢ, i★, j★, :), zc; label)
    scatterlines!(axN, interior(N², i★, j★, :), zf; label)
end

zm = -500

xlims!(axT, -2, 30)
xlims!(axS, 32, 40)
ylims!(axT, zm, 30)
ylims!(axS, zm, 30)
axislegend(axT, position=:rb)

display(fig)

