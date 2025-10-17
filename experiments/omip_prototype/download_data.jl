using ClimaOcean
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling: download_dataset

atmosphere = JRA55PrescribedAtmosphere(; dir="forcing_data/", 
                                         dataset=JRA55MultipleYears(), 
                                         backend=JRA55NetCDFBackend(10), 
                                         include_rivers_and_icebergs=true)