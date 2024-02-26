using NCDatasets
using GLMakie

filename = "single_column_omip_ocean_station_papa.nc"

ds = Dataset(filename)
K = ds["κᶜ"][1, 1, :, :]
zf = ds["zF"][:]
t = ds["time"][:]

# close(ds)

heatmap(t, zf, K')

