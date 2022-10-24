using ClimaOcean
using ClimaOcean.InterpolateData

degree   = 1/4
latitude = 75

bathymetry = interpolate_bathymetry_from_file(degree, latitude; passes = 20, interpolation_method = LinearInterpolation())
bathymetry = remove_connected_regions(bathymetry)

write_bathymetry_to_file!("bathymetry", bathymetry, latitude)