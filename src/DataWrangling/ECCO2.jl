module ECCO2

export ecco2_field, ecco2_center_mask, adjusted_ecco_tracers, initialize!

using ClimaOcean.InitialConditions: adjust_tracer!, three_dimensional_regrid!

using Oceananigans
using Oceananigans: architecture
using Oceananigans.BoundaryConditions
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index
using NCDatasets

temperature_filename = "THETA.1440x720x50.19920102.nc"
salinity_filename = "SALT.1440x720x50.19920102.nc"
effective_ice_thickness_filename = "SIheff.1440x720.19920102.nc"

ecco2_tracer_fields = Dict(
    :ecco2_temperature => :temperature,
    :ecco2_salinity => :salinity,
    :ecco_2_effective_ice_thickness => :effective_ice_thickness
)

ecco2_short_names = Dict(
    :temperature   => "THETA",
    :salinity      => "SALT",
    :effective_ice_thickness => "SIheff"
)

ecco2_location = Dict(
    :temperature   => (Center, Center, Center),
    :salinity      => (Center, Center, Center),
    :effective_ice_thickness => (Center, Center, Nothing)
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

    FT    = eltype(grid)
    field = Field{Center, Center, LZ}(grid)
    data  = convert.(FT, data)

    set!(field, data)
    fill_halo_regions!(field)

    return field
end

@kernel function _set_ecco2_mask!(mask, Tᵢ, minimum_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = ifelse(Tᵢ[i, j, k] < minimum_value, 0, 1)
end

function ecco2_center_mask(architecture = CPU(); minimum_value = Float32(-1e5))
    Tᵢ   = ecco2_field(:temperature; architecture)
    mask = CenterField(Tᵢ.grid)

    # Set the mask with ones where T is defined
    launch!(architecture, Tᵢ.grid, :xyz, _set_ecco2_mask!, mask, Tᵢ, minimum_value)

    return mask
end

function adjusted_ecco_field(variable_name; 
                             architecture = CPU(),
                             overwrite_existing = false, 
                             filename = "./data/initial_ecco_tracers.nc",
                             mask = ecco2_center_mask(architecture))
    
    if overwrite_existing || !isfile(filename)
        f = ecco2_field(variable_name; architecture)
        
        # Make sure all values are extended properly
        adjust_tracer!(f; mask)

        ds = Dataset(filename, "c")
        defVar(ds, string(variable_name), f, ("lat", "lon", "z"))

        close(ds)
    else
        ds = Dataset(filename)

        if haskey(ds, string(variable_name))
            f = ds[variable_name][:, :, :]
        else
            f = ecco2_field(variable_name; architecture)
            # Make sure all values are extended properly
            adjust_tracer!(f; mask)

            defVar(ds, string(variable_name), f, ("lat", "lon", "z"))
        end

        close(ds)
    end

    return f
end

function initialize!(model;
                     overwrite_existing = false,
                     filename = "./data/initial_ecco_tracers.nc", 
                     kwargs...)
    
    grid = model.grid
    arch = architecture(grid)
    
    custom_fields = Dict()
    ecco2_fields  = Dict()

    # Differentiate between custom initialization and ECCO2 initialization
    for (fldname, value) in kwargs
        if value ∈ keys(ecco2_tracer_fields)
            ecco2_fields[fldname] = ecco2_tracer_fields[value]
        else
            custom_fields[fldname] = value
        end
    end

    # Additional tracers not present in the ECCO dataset
    if !isempty(custom_fields)
        set!(model; custom_fields...)
    end

    # Set tracers from ECCO2
    if !isempty(ecco2_fields)
        mask = ecco2_center_mask(architecture(grid))
    
        for (fldname, variable_name) in ecco2_fields
            f = adjusted_ecco_field(variable_name; 
                                    architecture = arch,
                                    overwrite_existing, 
                                    filename,
                                    mask)

            f_grid = Field(ecco2_location[variable_name], grid)     
            three_dimensional_regrid!(f_grid, f)

            set!(model; fldname => f_grid)
        end
    end

    return nothing
end

end # module

