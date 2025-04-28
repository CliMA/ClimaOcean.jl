using ClimaOcean
using PythonCall

dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSDaily()
u_meta = ClimaOcean.DataWrangling.Metadatum(:u_velocity; dataset)
v_meta = ClimaOcean.DataWrangling.Metadatum(:v_velocity; dataset)
T_meta = ClimaOcean.DataWrangling.Metadatum(:temperature; dataset)
S_meta = ClimaOcean.DataWrangling.Metadatum(:salinity; dataset)

ClimaOcean.DataWrangling.download_dataset(u_meta)
ClimaOcean.DataWrangling.download_dataset(v_meta)
ClimaOcean.DataWrangling.download_dataset(T_meta)
ClimaOcean.DataWrangling.download_dataset(S_meta)
