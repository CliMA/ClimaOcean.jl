using ClimaOcean
using ClimaOcean.JRA55
using ClimaOcean.DataWrangling: download_dataset

atmosphere = JRA55PrescribedAtmosphere(; dir="forcing_data/", 
                                         dataset=MultiYearJRA55(), 
                                         backend=JRA55NetCDFBackend(100), 
                                         include_rivers_and_icebergs=true)
