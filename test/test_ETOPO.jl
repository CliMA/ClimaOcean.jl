using ClimaOcean.DataWrangling: download_dataset
using ClimaOcean.DataWrangling.ETOPO
using ClimaOcean

ETOPOmetadata = Metadatum(:temperature, dataset=ETOPOBathymetry())

download_dataset(ETOPOmetadata)