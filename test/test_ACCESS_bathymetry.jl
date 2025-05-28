using Oceananigans
using Statistics
using ClimaOcean

using ClimaOcean.Bathymetry: remove_minor_basins!
using ClimaOcean.DataWrangling.ACCESS

ACCESS_1dg_metadata = Metadatum(:bottom_height, dataset=ACCESS1deg())
ACCESS_025dg_metadata = Metadatum(:bottom_height, dataset=ACCESS025deg())

@show ACCESS_1dg_metadata
@show ACCESS_025dg_metadata