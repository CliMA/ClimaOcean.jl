using ClimaOcean.DataWrangling: download_dataset
using ClimaOcean.DataWrangling.ETOPO
using ClimaOcean
using NCDatasets
using Oceananigans


data_path = expanduser("/Users/tsohail/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/uom/ocean-ensembles/data/")

ETOPOmetadata = Metadatum(:bottom_height, dataset=ETOPOBathymetry(), dir = data_path)

download_dataset(ETOPOmetadata)

filepath = "/Users/tsohail/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/uom/ocean-ensembles/data/ETOPO_2022_v1_60s_N90W180_surface.nc"
dataset = Dataset(filepath, "r")

Nx = Integer(360)
Ny = Integer(180)
Nz = Integer(50)

arch = CPU()

z_faces = (-4000, 0)

underlying_grid = TripolarGrid(arch;
                               size = (Nx, Ny, Nz),
                               z = z_faces,
                               halo = (5, 5, 4),
                               first_pole_longitude = 70,
                               north_poles_latitude = 55)

bottom_height = regrid_bathymetry(underlying_grid, ETOPOmetadata;
                                  minimum_depth = 10,
                                  interpolation_passes = 75, # 75 interpolation passes smooth the bathymetry near Florida so that the Gulf Stream is able to flow
				                  major_basins = 2)                              