using ClimaOcean
using ClimaOcean.Diagnostics: MixedLayerDepthField
using ClimaOcean.DataWrangling.ECCO: ECCOFieldTimeSeries
using Oceananigans
using CairoMakie
using Printf
using Dates

using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans.BuoyancyFormulations: buoyancy

arch = CPU()
Nx = 360
Ny = 160

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

start_date = DateTime(1993, 1, 1)
end_date   = DateTime(2003, 1, 1)

Tt = ECCOFieldTimeSeries(:temprature, grid; start_date, end_date, time_indices_in_memory=2)
St = ECCOFieldTimeSeries(:salinity,   grid; start_date, end_date, time_indices_in_memory=2)
ht = FieldTimeSeries{Center, Center, Nothing}(grid, Tt.times)

equation_of_state = TEOS10EquationOfState()
sb = SeawaterBuoyancy(; equation_of_state)
tracers = (T=Tt[1], S=St[1])
h = MixedLayerDepthField(sb, grid, tracers)

Nt = length(ht)
for n = 1:Nt-1
    local tracers
    tracers = (T=Tt[n], S=St[n])
    h.operand.buoyancy_perturbation = buoyancy(sb, grid, tracers)
    @show n
    @time compute!(h)
    parent(ht[n]) .= parent(h)
end

function titlestr(n)
    d = dates[n]
    yr = year(d)
    mn = monthname(d)
    return string("ECCO mixed layer depth on ", mn, " ", yr)
end

fig = Figure(size=(1500, 800))
axh = Axis(fig[2, 1], xlabel="Longitude", ylabel="Latitude")
n = Observable(1)

str = @lift titlestr($n)
Label(fig[1, 1], str, tellwidth=false)

hn = @lift ht[$n]
hm = heatmap!(axh, hn, colorrange=(0, 500), colormap=:magma, nan_color=:lightgray)
Colorbar(fig[2, 2], hm, label="Mixed layer depth (m)")
fig

# And record a movie

record(fig, "ecco_mld.mp4", 1:Nt-1, framerate=4) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
