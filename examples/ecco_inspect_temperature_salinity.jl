using Oceananigans
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

using CairoMakie
using Printf
using ClimaOcean
using ClimaOcean.DataWrangling.ECCO: ECCO_field, ECCOFieldTimeSeries, ECCO4Monthly
using CFTime
using Dates

arch = CPU()
Nx = 360 ÷ 4
Ny = 160 ÷ 4

z = ClimaOcean.DataWrangling.ECCO.ECCO_z
z = z[20:end]
Nz = length(z) - 1

grid = LatitudeLongitudeGrid(arch; z,
                             size = (Nx, Ny, Nz),
                             latitude  = (-80, 80),
                             longitude = (0, 360))

bottom_height = regrid_bathymetry(grid; 
                                  minimum_depth = 10,
                                  interpolation_passes = 5,
                                  major_basins = 1)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

T = CenterField(grid)
S = CenterField(grid)

using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans.BuoyancyFormulations: buoyancy

equation_of_state = TEOS10EquationOfState()
sb = SeawaterBuoyancy(; equation_of_state)
tracers = (T=T, S=S)
b = Field(buoyancy(sb, grid, tracers))

start = DateTime(1993, 1, 1)
stop  = DateTime(1999, 1, 1)
dates = range(start; stop, step=Month(1))

Tmeta = Metadata(:temperature; dates, dataset=ECCO4Monthly())
Smeta = Metadata(:salinity; dates, dataset=ECCO4Monthly())

Tt = ECCOFieldTimeSeries(Tmeta, grid; time_indices_in_memory=length(dates))
St = ECCOFieldTimeSeries(Smeta, grid; time_indices_in_memory=length(dates))

fig = Figure(size=(900, 1050))

axT = Axis(fig[1, 1])
axS = Axis(fig[2, 1])
axb = Axis(fig[3, 1])

Nt = length(dates)
grid = T.grid
Nz = size(grid, 3)
kslider = Slider(fig[1:3, 0], range=1:Nz, startvalue=Nz, horizontal=false)
nslider = Slider(fig[4, 1:2], range=1:Nt, startvalue=1)
k = kslider.value
n = nslider.value

Tk = @lift view(Tt[$n], :, :, $k)
Sk = @lift view(St[$n], :, :, $k)

Δb = @lift begin
    parent(T) .= parent(Tt[$n])
    parent(S) .= parent(St[$n])
    compute!(b)
    mask_immersed_field!(b, NaN)
    Δb = interior(b, :, :, Nz) .- interior(b, :, :, $k)
    Δb
end

hmT = heatmap!(axT, Tk, nan_color=:lightgray, colorrange=(-2, 30), colormap=:thermal)
hmS = heatmap!(axS, Sk, nan_color=:lightgray, colorrange=(31, 37), colormap=:haline)
hmb = heatmap!(axb, Δb, nan_color=:lightgray, colorrange=(0, 1e-3), colormap=:magma)

Colorbar(fig[1, 2], hmT, label="Temperature (ᵒC)")
Colorbar(fig[2, 2], hmS, label="Salinity (g kg⁻¹)")
Colorbar(fig[3, 2], hmb, label="Buoyancy difference (m s⁻²)")

fig
