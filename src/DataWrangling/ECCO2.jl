module ECCO2

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
    :salinity      => "DEPTH_S",
)

variable_is_three_dimensional = Dict(
    :temperature             => true,
    :salinity                => true,
    :effective_ice_thickness => false,
)

ecco2_file_names = Dict(
    :temperature             => "ecco2_temperature_19920112.nc",
    :salinity                => "ecco2_salinity_19920112.nc",
    :effective_ice_thickness => "ecco2_effective_ice_thickness_19920112.nc",
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
                     halo = (1, 1, 1),
                     url = ecco2_urls[variable_name],
                     filename = ecco2_file_names[variable_name],
                     short_name = ecco2_short_names[variable_name])

    isfile(filename) || download(url, filename)

    ds = Dataset(filename)

    longitude = (0, 360)
    latitude = (-90, 90)
    topology = (Periodic, Bounded, Bounded)

    if variable_is_three_dimensional[variable_name] 
        data = ds[short_name][:, :, :, 1]
        depth_name = ecco2_depth_names[variable_name]
        z = construct_vertical_interfaces(ds, depth_name)
        Nx, Ny, Nz = size(data)
        LZ = Center
    else
        data = ds[short_name][:, :, 1]
        Nx, Ny = size(data)
        z = nothing
        LZ = Nothing
    end

    close(ds)

    grid = LatitudeLongitudeGrid(architecture; halo, size, topology,
                                 longitude, latitude, z)

    field = Field{Center, Center, LZ}(grid)
    set!(field, data)
    fill_halo_regions!(field)

    return field
end

function ecco2_bottom_height_from_temperature()
    Tᵢ = ecco2_field(:temperature, CPU())

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

