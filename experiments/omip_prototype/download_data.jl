using ClimaOcean
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling: download_dataset

dir = "forcing_data/"

atmosphere = JRA55PrescribedAtmosphere(; dataset=JRA55MultipleYears(), dir, include_rivers_and_icebergs=true)