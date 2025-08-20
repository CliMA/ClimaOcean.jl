using ClimaOcean
using PythonCall
using Oceananigans

VerosModule = Base.get_extension(ClimaOcean, :ClimaOceanPythonCallExt)
VerosModule.remove_outputs(:global_4deg)

ocean = VerosModule.veros_ocean_simulation("global_4deg", :GlobalFourDegreeSetup)
VerosModule.veros_settings_set!(ocean, "dt_tracer", 1800.0)

atmos = JRA55PrescribedAtmosphere(; backend = JRA55NetCDFBackend(10))
radiation = Radiation()
coupled_model = OceanSeaIceModel(ocean, nothing; atmosphere=atmos, radiation)

