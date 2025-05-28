include("runtests_setup.jl")

using Oceananigans
using Statistics
using ClimaOcean

using ClimaOcean.Bathymetry: remove_minor_basins!
using ClimaOcean.DataWrangling.ACCESS
# using ClimaOcean.DataWrangling.ETOPO

# ETOPO_meta = Metadatum(:bottom_height, dataset=ETOPO2022())
data_path = expanduser("/Users/tsohail/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/uom/ocean-ensembles/data/")

ACCESS_1dg_metadata = Metadatum(:bottom_height, dataset=ACCESS1deg(), dir = data_path)
ACCESS_025dg_metadata = Metadatum(:bottom_height, dataset=ACCESS025deg(), dir = data_path)

@show ACCESS_1dg_metadata
@show ACCESS_025dg_metadata

# Testing downloading
ClimaOcean.DataWrangling.download_dataset(ACCESS_1dg_metadata)
@test isfile(metadata_path(ACCESS_1dg_metadata))
# Testing downloading
ClimaOcean.DataWrangling.download_dataset(ACCESS_025dg_metadata)
@test isfile(metadata_path(ACCESS_025dg_metadata))

