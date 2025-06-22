using ClimaOcean
using Oceananigans
using PythonCall

arch = CPU()
Nx = 20 * 12
Ny = 20 * 12
Nz = 50

depth = 6000
zgrid = exponential_vertical_faces(; Nz, depth, scale=depth/4.5)
zf = z_faces(zgrid)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces,
                             latitude  = (35, 55),
                             longitude = (200, 220))

bounding_box = ClimaOcean.DataWrangling.BoundingBox(longitude=(200, 220), latitude=(35, 55))

# dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSStatic()
# static_meta = ClimaOcean.DataWrangling.Metadatum(:depth; dataset, bounding_box)
# coords_path = ClimaOcean.DataWrangling.download_dataset(static_meta)
# @info "Downloaded coordinates data to $coords_path"

# T_ecco = ClimaOcean.DataWrangling.ECCOMetadatum(:temperature; dataset, bounding_box)
# T_en4_meta = ClimaOcean.DataWrangling.EN4Metadatum(:temperature)
# T_en4_path = ClimaOcean.DataWrangling.download_dataset(T_en4_meta)
# T_en4 = Field(T_en4_meta)

dataset = ClimaOcean.DataWrangling.Copernicus.GLORYSDaily()
T_meta = ClimaOcean.DataWrangling.Metadatum(:temperature; dataset, bounding_box)
T_path = ClimaOcean.DataWrangling.download_dataset(T_meta)
@info "Downloaded temperature data to $T_path"
T = Field(T_meta, inpainting=nothing)

#=
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
=#

#=
# FOR ERA5:
# account: https://cds.climate.copernicus.eu/how-to-api
CondaPkg.add("cdsapi"; channel = "conda-forge")
cds = pyimport("cdsapi")          # should succeed instantly
client = cds.Client()
=#
