module ECCO2

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!

using NCDatasets

const ECCO2_Nx = 1440
const ECCO2_Ny = 720
const ECCO2_Nz = 50

# Vertical coordinate
const ECCO2_z = [
    -6128.75,
    -5683.75,
    -5250.25,
    -4839.75,
    -4452.25,
    -4087.75,
    -3746.25,
    -3427.75,
    -3132.25,
    -2859.75,
    -2610.25,
    -2383.74,
    -2180.13,
    -1999.09,
    -1839.64,
    -1699.66,
    -1575.64,
    -1463.12,
    -1357.68,
    -1255.87,
    -1155.72,
    -1056.53,
    -958.45,
    -862.10,
    -768.43,
    -678.57,
    -593.72,
    -515.09,
    -443.70,
    -380.30,
    -325.30,
    -278.70,
    -240.09,
    -208.72,
    -183.57,
    -163.43,
    -147.11,
    -133.45,
    -121.51,
    -110.59,
    -100.20,
    -90.06,
    -80.01,
    -70.0,
    -60.0,
    -50.0,
    -40.0,
    -30.0,
    -20.0,
    -10.0,
      0.0,
]

filenames = Dict(
    :temperature       => "THETA.1440x720x50.19920102.nc",
    :salinity          => "SALT.1440x720x50.19920102.nc",
    :sea_ice_thickness => "SIheff.1440x720.19920102.nc",
)

shortnames = Dict(
    :temperature       => "THETA",
    :salinity          => "SALT",
    :sea_ice_thickness => "SIheff",
)

depthnames = Dict(
    :temperature       => "DEPTH_T",
    :salinity          => "DEPTH_S",
    :sea_ice_thickness => nothing,
)

# Downloaded from https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N
# These files are just for Jan 2, 1992.
urls = Dict(
    :temperature => "https://www.dropbox.com/scl/fi/01h96yo2fhnnvt2zkmu0d/" *
                    "THETA.1440x720x50.19920102.nc?rlkey=ycso2v09gc6v2qb5j0lff0tjs&dl=0",

    :salinity => "https://www.dropbox.com/scl/fi/t068we10j5skphd461zg8/" *
                 "SALT.1440x720x50.19920102.nc?rlkey=r5each0ytdtzh5icedvzpe7bw&dl=0",

    :sea_ice_thickness => "https://www.dropbox.com/scl/fi/x0v9gjrfebwsef4tv1dvn/" *
                          "SIheff.1440x720.19920102.nc?rlkey=2uel3jtzbsplr28ejcnx3u6am&dl=0",
)

surface_variable(variable_name) = variable_name == :sea_ice_thickness

function ecco2_field(variable_name;
                     architecture = CPU(),
                     filename = filenames[variable_name],
                     shortname = shortnames[variable_name],
                     depthname = depthnames[variable_name],
                     url = urls[variable_name])
                     
    isfile(filename) || download(url, filename)

    ds = Dataset(filename)

    grid = LatitudeLongitudeGrid(architecture,
                                 size = (ECCO2_Nx, ECCO2_Ny, ECCO2_Nz),
                                 longitude = (0, 360),
                                 latitude = (-90, 90),
                                 z = ECCO2_z,
                                 halo = (1, 1, 1),
                                 topology = (Periodic, Bounded, Bounded))

    if surface_variable(variable_name)
        field = Field{Center, Center, Nothing}(grid)
        data = ds[shortname][:, :, 1]
        data = convert(Array{Float32, 2}, data)
    else
        field = CenterField(grid) # u, v not supported
        data = ds[shortname][:, :, :, 1]
        data = convert(Array{Float32, 3}, data)
        data = reverse(data, dims=3)
    end

    set!(field, data)
    fill_halo_regions!(field)

    return field
end

end # module
