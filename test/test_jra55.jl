using Oceananigans
using Oceananigans.Fields: xnode, ynode, znode, interpolate
using NCDatasets
using Downloads: download

shortwave_radiation_filename = "RYF.rsds.1990_1991.nc"
shortwave_radiation_url = "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                          "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0"

isfile(shortwave_radiation_filename) || download(shortwave_radiation_url, shortwave_radiation_filename)

shortwave_radiation_ds = Dataset(shortwave_radiation_filename)
shortwave_radiation = shortwave_radiation_ds["rsds"][:, :, :]

# Make source field
λ = shortwave_radiation_ds["lon_bnds"][1, :]
φ = shortwave_radiation_ds["lat_bnds"][1, :]
close(shortwave_radiation_ds)

push!(φ, 90)
push!(λ, λ[1] + 360)

Nxs = length(λ) - 1
Nys = length(φ) - 1

source_grid = LatitudeLongitudeGrid(size = (Nxs, Nys, 1);
                                    longitude = λ,
                                    latitude = φ,
                                    z = (0, 1),
                                    topology = (Periodic, Bounded, Bounded))

source_field = Field{Center, Center, Center}(source_grid)
set!(source_field, shortwave_radiation[:, :, 1])

# Make target grid and field
# Quarter degree target grid with some latitude bounds
Nxt = 1440

southern_limit = -79
northern_limit = -30
j₁ = 4 * (90 + southern_limit)
j₂ = 720 - 4 * (90 - northern_limit) + 1
Nyt = j₂ - j₁ + 1

target_grid = LatitudeLongitudeGrid(size = (Nxt, Nyt, 1);
                                    longitude = (0, 360),
                                    latitude = (southern_limit, northern_limit),
                                    z = (0, 1),
                                    topology = (Periodic, Bounded, Bounded))

Qˢʷ = target_field = Field{Center, Center, Center}(target_grid)

# Make target grid and field
source_location = Tuple(L() for L in location(source_field))
target_location = Tuple(L() for L in location(target_field))

for i = 1:Nxt, j = 1:Nyt
    k = 1
    x = xnode(i, j, k, target_grid, target_location...)
    y = ynode(i, j, k, target_grid, target_location...)
    z = 0.5 #znode(i, j, k, grid, target_location...)
    @inbounds target_field[i, j, k] = interpolate(source_field, source_location..., source_grid, x, y, z)
end

source_grid = RectilinearGrid(size=(2, 2), x=(0, 1), y=(0, 1), topology=(Periodic, Periodic, Flat))

s = Field{Center, Center, Nothing}(source_grid)
set!(s, (x, y) -> x + y)

loc = Tuple(L() for L in location(s))

x = y = 0.67
z = 0

interpolate(s, loc..., source_grid, x, y, z)

