module ECCO2

using Dates
using Oceananigans
using Oceananigans.BoundaryConditions
using NCDatasets

temperature_filename = "THETA.1440x720x50.19920102.nc"
salinity_filename = "SALT.1440x720x50.19920102.nc"
effective_ice_thickness_filename = "SIheff.1440x720.19920102.nc"

ecco2_short_names = Dict(
    :temperature   => "THETA",
    :salinity      => "SALT",
    :effective_ice_thickness => "SIheff"
)

ecco2_depth_names = Dict(
    :temperature   => "DEPTH_T",
    :salinity      => "DEPTH_T",
)

variable_is_three_dimensional = Dict(
    :temperature             => true,
    :salinity                => true,
    :effective_ice_thickness => false,
)

ecco2_file_names = Dict(
    :temperature             => "ecco2_temperature_19920102.nc",
    :salinity                => "ecco2_salinity_19920102.nc",
    :effective_ice_thickness => "ecco2_effective_ice_thickness_19920102.nc",
)

# Downloaded from https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N

ecco2_urls = Dict(
    :temperature => "https://www.dropbox.com/scl/fi/01h96yo2fhnnvt2zkmu0d/" *
                    "THETA.1440x720x50.19920102.nc?rlkey=ycso2v09gc6v2qb5j0lff0tjs&dl=0",

    :salinity => "https://www.dropbox.com/scl/fi/t068we10j5skphd461zg8/" *
                 "SALT.1440x720x50.19920102.nc?rlkey=r5each0ytdtzh5icedvzpe7bw&dl=0",

    :effective_ice_thickness => "https://www.dropbox.com/scl/fi/x0v9gjrfebwsef4tv1dvn/" *
                                "SIheff.1440x720.19920102.nc?rlkey=2uel3jtzbsplr28ejcnx3u6am&dl=0"
)

function construct_vertical_interfaces(ds, depth_name)
    # Construct vertical coordinate
    depth = ds[depth_name][:]
    zc = -reverse(depth)

    # Interface depths from cell center depths
    zf = (zc[1:end-1] .+ zc[2:end]) ./ 2
    push!(zf, 0)
    
    Δz = zc[2] - zc[1]
    pushfirst!(zf, zf[1] - Δz)

    return zf
end

function ecco2_field(variable_name;
                     architecture = CPU(),
                     horizontal_halo = (1, 1),
                     url = ecco2_urls[variable_name],
                     filename = ecco2_file_names[variable_name],
                     short_name = ecco2_short_names[variable_name])

    isfile(filename) || download(url, filename)

    ds = Dataset(filename)

    longitude = (0, 360)
    latitude = (-90, 90)
    TX, TY = (Periodic, Bounded)

    if variable_is_three_dimensional[variable_name] 
        data = ds[short_name][:, :, :, 1]
        depth_name = ecco2_depth_names[variable_name]
        
        # The surface layer in three-dimensional ECCO fields is at `k = 1`
        data = reverse(data, dims = 3)
        
        z    = construct_vertical_interfaces(ds, depth_name)
        N    = size(data)

        # add vertical halo for 3D fields
        halo = (horizontal_halo..., 1)

        LZ   = Center
        TZ   = Bounded
    else
        data = ds[short_name][:, :, 1]
        N    = size(data)
        z    = nothing
        halo = horizontal_halo
        LZ   = Nothing
        TZ   = Flat
    end

    close(ds)

    # Flat in z if the variable is two-dimensional
    grid = LatitudeLongitudeGrid(architecture; halo, size = N, topology = (TX, TY, TZ),
                                 longitude, latitude, z)

    field = Field{Center, Center, LZ}(grid)
    
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

# Compute the dates that are available
first_date = Date(1992, 1, 2)
last_date = Date(2023, 4, 28)
data_frequency = 3
days_available = Day(last_date - first_date)
Ndates = Int(days_available.value / 3)
available_dates = [first_date + Day(3n) for n = 0:Ndates]

function wget_ecco2_data(variable_name, user_name, password;
                         url = "https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N"
                         date = first_date)

    # Try to be helpful
    warnmsg = string(date,
                     " may not be available! As far as we know, data is available", '\n',
                     " every 3 days starting from 1992-01-02".)


    date ∈ available_dates || @warn(warnmsg)

    date_str = Dates.format(date, "yyyymmdd")
    field_name = ecco2_short_names[variable_name]
    filename = string(field_name, ".1440x720.", date_str, ".nc")

    dirname = string(field_name, ".nc")
    fileurl = joinpath(url, dirname, filename)
    cmd = `wget --http-user=$(user_name) --http-passwd=$(password) $(fileurl)`
    run(cmd)

    return filename
end

end # module

