module ECCO2

using Dates

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!

using NCDatasets

# Data from
#
# https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N
#
# These files are just for Jan 2, 1992.

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

filenames_19920102 = Dict(
    :temperature           => "THETA.1440x720x50.19920102.nc",
    :salinity              => "SALT.1440x720x50.19920102.nc",
    :sea_ice_thickness     => "SIheff.1440x720.19920102.nc",
    :sea_ice_area_fraction => "SIarea.1440x720.19920102.nc",
    :u_velocity            => "UVEL.1440x720.19920102.nc",
    :v_velocity            => "VVEL.1440x720.19920102.nc",
)

filenames_19921001 = Dict(
    :temperature           => "THETA.1440x720x50.19921001.nc",
    :salinity              => "SALT.1440x720x50.19921001.nc",
    :sea_ice_thickness     => "SIheff.1440x720.19921001.nc",
    :sea_ice_area_fraction => "SIarea.1440x720.19921001.nc",
    :u_velocity            => "UVEL.1440x720.19921001.nc",
    :v_velocity            => "VVEL.1440x720.19921001.nc",
)

urls_19920102 = Dict(
    :temperature           => "https://www.dropbox.com/scl/fi/01h96yo2fhnnvt2zkmu0d/THETA.1440x720x50.19920102.nc?rlkey=ycso2v09gc6v2qb5j0lff0tjs",
    :salinity              => "https://www.dropbox.com/scl/fi/t068we10j5skphd461zg8/SALT.1440x720x50.19920102.nc?rlkey=r5each0ytdtzh5icedvzpe7bw",
    :sea_ice_thickness     => "https://www.dropbox.com/scl/fi/x0v9gjrfebwsef4tv1dvn/SIheff.1440x720.19920102.nc?rlkey=2uel3jtzbsplr28ejcnx3u6am",
    :sea_ice_area_fraction => "https://www.dropbox.com/scl/fi/q14moq3201zicppu8ff8h/SIarea.1440x720.19920102.nc?rlkey=pt7pt80gr7r6mmjm9e0u4f5n1",
    :u_velocity            => "https://www.dropbox.com/scl/fi/myur9kpanc5mprrf5ge32/UVEL.1440x720x50.19920102.nc?rlkey=7a5dpvfgoc87yr6q5ktrqwndu",
    :v_velocity            => "https://www.dropbox.com/scl/fi/buic35gssyeyfqohenkeo/VVEL.1440x720x50.19920102.nc?rlkey=fau48w4t5ruop4s6gm8t7z0a0",
)

urls_19921001 = Dict(
    :temperature           => "https://www.dropbox.com/scl/fi/169f3981460uhk9h69k0f/THETA.1440x720x50.19921001.nc?rlkey=mgal3xt0qy2c59y395ybio11v",
    :salinity              => "https://www.dropbox.com/scl/fi/f9zfm34vqz732jrrhjrg3/SALT.1440x720x50.19921001.nc?rlkey=y5dv0s41gb6f9guvu0iorw28p",
    :sea_ice_thickness     => "https://www.dropbox.com/scl/fi/mtmziurepom8kpjn82d07/SIheff.1440x720.19921001.nc?rlkey=9uhuxg2n9iw6894afj4t53drv",
    :sea_ice_area_fraction => "https://www.dropbox.com/scl/fi/ntflhyrmsnit9vco402co/SIarea.1440x720.19921001.nc?rlkey=eakzc788btql1q6ndj9l8cr2q",
    #:u_velocity            => "https://www.dropbox.com/scl/fi/e6s9c013r2ddift4f8ugi/UVEL.1440x720x50.19921001.nc?rlkey=fpd7mv1zv3fkmyg8w11b94sbp&dl=0",
    :u_velocity            => "https://www.dropbox.com/scl/fi/e6s9c013r2ddift4f8ugi/UVEL.1440x720x50.19921001.nc?rlkey=fpd7mv1zv3fkmyg8w11b94sbp&dl=0",
    :v_velocity            => "https://www.dropbox.com/scl/fi/nxuohvhvdu0ig552osf1d/VVEL.1440x720x50.19921001.nc?rlkey=vz4ttp3myxhertdxvt1lyjp1d",
)

filenames = Dict(
    "1992-01-02" => filenames_19920102,
    "1992-10-01" => filenames_19921001,
)

urls = Dict(
    "1992-01-02" => urls_19920102,
    "1992-10-01" => urls_19921001,
)

shortnames = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_area_fraction => "SIarea",
    :u_velocity            => "UVEL",
    :v_velocity            => "VVEL",
)



surface_variable(variable_name) = variable_name == :sea_ice_thickness

function ecco2_field(variable_name, date=Date(1992, 01, 02);
                     architecture = CPU(),
                     filename  = filenames[string(date)][variable_name],
                     url       =      urls[string(date)][variable_name],
                     shortname = shortnames[variable_name])
                     
    if !isfile(filename)
        print("Downloading $filename...")
        start_time = time_ns()
        download(url, filename)
        elapsed = time_ns() - start_time
        print(" done (", prettytime(elapsed * 1e-9), ").", '\n')
    end

    ds = Dataset(filename)

    grid = LatitudeLongitudeGrid(architecture,
                                 size = (ECCO2_Nx, ECCO2_Ny, ECCO2_Nz),
                                 longitude = (0, 360),
                                 latitude = (-90, 90),
                                 z = ECCO2_z,
                                 halo = (1, 1, 1),
                                 topology = (Periodic, Bounded, Bounded))

    # TODO: figure out what's going on with the locations
    if surface_variable(variable_name)
        field = Field{Center, Center, Nothing}(grid)
        data = ds[shortname][:, :, 1]
        data = convert(Array{Float32, 2}, data)
    else
        field = CenterField(grid)
        data = ds[shortname][:, :, :, 1]
        data = convert(Array{Float32, 3}, data)
        data = reverse(data, dims=3)
    end

    set!(field, data)
    fill_halo_regions!(field)

    return field
end

function ecco2_bottom_height_from_temperature()
    Tᵢ = ecco2_field(:temperature)

    missing_value = Float32(-9.9e22)

    # Construct bottom_height depth by analyzing T
    Nx, Ny, Nz = size(Tᵢ)
    bottom_height = ones(Nx, Ny) .* (zf[1] - Δz)
    zf = znodes(Tᵢ.grid, Face())
    
    for i = 1:Nx, j = 1:Ny
        @inbounds for k = Nz:-1:1
            if Tᵢ[i, j, k] < -10
                bottom_height[i, j] = zf[k+1]
                break
            end
        end
    end

    return bottom_height
end

end # module

