using ClimaOcean
using NCDatasets

cachepath = ClimaOcean.DataWrangling.JRA55.download_jra55_cache
filename = "RYF.tas.1990_1991.nc"
filepath = joinpath(cachepath, filename)

ds = Dataset(filepath)
Nx, Ny, Nt = size(ds["tas"])
ds["tas"][1, 1, [Nt, 1]]
close(ds)
