using ClimaOcean
using Oceananigans
using PythonCall

dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSStatic()
static_meta = ClimaOcean.DataWrangling.Metadatum(:depth; dataset)
coords_path = ClimaOcean.DataWrangling.download_dataset(static_meta)
@info "Downloaded coordinates data to $coords_path"

dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSDaily()
T_meta = ClimaOcean.DataWrangling.Metadatum(:temperature; dataset)
T_path = ClimaOcean.DataWrangling.download_dataset(T_meta)
@info "Downloaded temperature data to $T_path"
T = Field(T_meta)

u_meta = ClimaOcean.DataWrangling.Metadatum(:u_velocity; dataset)
u_path = ClimaOcean.DataWrangling.download_dataset(u_meta)
@info "Downloaded u velocity data to $u_path"
u = Field(u_meta)

v_meta = ClimaOcean.DataWrangling.Metadatum(:v_velocity; dataset)
v_path = ClimaOcean.DataWrangling.download_dataset(v_meta)
@info "Downloaded data to $v_path"
v = Field(v_meta)

S_meta = ClimaOcean.DataWrangling.Metadatum(:salinity; dataset)
S_path = ClimaOcean.DataWrangling.download_dataset(S_meta)
@info "Downloaded data to $S_path"
S = Field(S_meta)