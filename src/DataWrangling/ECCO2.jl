module ECCO2

export ECCO2Metadata, ecco2_field, ecco2_center_mask, adjusted_ecco_tracers, initialize!

using ClimaOcean.DataWrangling: inpaint_mask!
using ClimaOcean.InitialConditions: three_dimensional_regrid!, three_dimensional_interpolate!

using Oceananigans
using Oceananigans.Architectures: architecture, child_architecture
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations: DistributedField, all_reduce, barrier!
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index
using NCDatasets
using Downloads: download

import Oceananigans.Fields: set!

# Ecco field used to set model's initial conditions
struct ECCO2Metadata
    name  :: Symbol
    year  :: Int
    month :: Int
    day   :: Int
end

# We only have 1992 at the moment
ECCO2Metadata(name::Symbol) = ECCO2Metadata(name, 1992, 1, 2)

filename(data::ECCO2Metadata) = "ecco2_" * string(data.name) * "_$(data.year)$(data.month)$(data.day).nc"

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

surface_variable(variable_name) = variable_name == :sea_ice_thickness

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

function empty_ecco2_field(data::ECCO2Metadata;
                           architecture = CPU(), 
                           horizontal_halo = (5, 5))

    variable_name = data.name

    location = ecco2_location[variable_name]

    longitude = (0, 360)
    latitude = (-90, 90)
    TX, TY = (Periodic, Bounded)

    filename = ecco2_file_names[variable_name]
    
    ds = Dataset(filename)

    if variable_is_three_dimensional[variable_name] 
        depth_name = ecco2_depth_names[variable_name]
        z    = construct_vertical_interfaces(ds, depth_name)
        # add vertical halo for 3D fields
        halo = (horizontal_halo..., 1)
        LZ   = Center
        TZ   = Bounded
        N    = (1440, 720, 50)
    else
        z    = nothing
        halo = horizontal_halo
        LZ   = Nothing
        TZ   = Flat
        N    = (1440, 720)
    end

    # Flat in z if the variable is two-dimensional
    grid = LatitudeLongitudeGrid(architecture; halo, size = N, topology = (TX, TY, TZ),
                                 longitude, latitude, z)

    return Field{location...}(grid)
end

"""
    ecco2_field(variable_name;
                architecture = CPU(),
                horizontal_halo = (1, 1),
                user_data = nothing,
                url = ecco2_urls[variable_name],
                filename = ecco2_file_names[variable_name],
                short_name = ecco2_short_names[variable_name])

Retrieve the ecco2 field corresponding to `variable_name`. 
The data is either:
(1) retrieved from `filename`,
(2) dowloaded from `url` if `filename` does not exists,
(3) filled from `user_data` if `user_data` is provided.
"""
function ecco2_field(variable_name;
                     architecture = CPU(),
                     horizontal_halo = (5, 5),
                     user_data = nothing,
                     year = 1992,
                     month = 1,
                     day  = 2, 
                     url  = ecco2_urls[variable_name],
                     filename = ecco2_file_names[variable_name],
                     short_name = ecco2_short_names[variable_name])

    ecco2_data = ECCO2Metadata(variable_name, year, month, day)

    isfile(filename) || download(url, filename)

    if user_data isa Nothing
        ds = Dataset(filename)
        
        if variable_is_three_dimensional[variable_name] 
            data = ds[short_name][:, :, :, 1]
            # The surface layer in three-dimensional ECCO fields is at `k = 1`
            data = reverse(data, dims = 3)
        else
            data = ds[short_name][:, :, 1]
        end        
    else
        data = user_data
    end

    field = empty_ecco2_field(ecco2_data; architecture, horizontal_halo)
    FT    = eltype(field)
    data  = convert.(FT, data)

    set!(field, data)
    fill_halo_regions!(field)

    return field
end

@kernel function _set_ecco2_mask!(mask, Tᵢ, minimum_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = Tᵢ[i, j, k] < minimum_value
end

"""
    ecco2_center_mask(architecture = CPU(); minimum_value = Float32(-1e5))

A boolean field where `false` represents a missing value in the ECCO2 :temperature dataset.
"""
function ecco2_center_mask(architecture = CPU(); minimum_value = Float32(-1e5))
    Tᵢ   = ecco2_field(:temperature; architecture)
    mask = CenterField(Tᵢ.grid, Bool)

    # Set the mask with ones where T is defined
    launch!(architecture, Tᵢ.grid, :xyz, _set_ecco2_mask!, mask, Tᵢ, minimum_value)

    return mask
end

"""
    inpainted_ecco2_field(variable_name; 
                          architecture = CPU(),
                          filename = "./inpainted_ecco2_fields.nc",
                          mask = ecco2_center_mask(architecture))
    
Retrieve the ECCO2 field corresponding to `variable_name` inpainted to fill all the
missing values in the original dataset.

Arguments:
==========

- `variable_name`: the variable name corresponding to the Dataset.

Keyword Arguments:
==================

- `architecture`: either `CPU()` or `GPU()`.

- `filename`: the path where to retrieve the data from. If the file does not exist,
              the data will be retrived from the ECCO2 dataset, inpainted, and
              saved to `filename`.

- `mask`: the mask used to inpaint the field (see `inpaint_mask!`).
"""
function inpainted_ecco2_field(variable_name; 
                               architecture = CPU(),
                               filename = "./inpainted_ecco2_fields.nc",
                               mask = ecco2_center_mask(architecture),
                               kw...)
    
    if !isfile(filename)
        f = ecco2_field(variable_name; architecture)
        
        # Make sure all values are extended properly
        @info "In-painting ecco $variable_name and saving it in $filename"
        inpaint_mask!(f, mask; kw...)

        fill_halo_regions!(f)

        ds = Dataset(filename, "c")
        defVar(ds, string(variable_name), Array(interior(f)), ("lat", "lon", "z"))

        close(ds)
    else
        ds = Dataset(filename, "a")

        if haskey(ds, string(variable_name))
            data = ds[variable_name][:, :, :]
            f = ecco2_field(variable_name; architecture, user_data = data)
            fill_halo_regions!(f)
        else
            f = ecco2_field(variable_name; architecture)
            # Make sure all values are inpainted properly
            @info "In-painting ecco $variable_name and saving it in $filename"
            inpaint_mask!(f, mask; kw...)
            fill_halo_regions!(f)

            defVar(ds, string(variable_name), Array(interior(f)), ("lat", "lon", "z"))
        end

        close(ds)
    end

    return f
end

function set!(field::DistributedField, ecco2_metadata::ECCO2Metadata; filename="./inpainted_ecco2_fields.nc", kw...)
    # Fields initialized from ECCO2
    grid = field.grid
    arch = architecture(grid)
    child_arch = child_architecture(arch)
    name = ecco2_metadata.name

    mask = ecco2_center_mask(child_arch)

    f_ecco = if arch.local_rank == 0 # Make sure we read/write the file using only one core
            inpainted_ecco2_field(name; filename, mask,
                                  architecture = child_arch,
                                  kw...)
    else
        empty_ecco2_field(ecco2_metadata; architecture = child_arch)
    end

    barrier!(arch)

    # Distribute ecco field to all workers
    parent(f_ecco) .= all_reduce(+, parent(f_ecco), arch)

    f_grid = Field(ecco2_location[name], grid)   
    three_dimensional_interpolate!(f_grid, f_ecco)
    set!(field, f_grid)
    
    return field
end

function set!(field::Field, ecco2_metadata::ECCO2Metadata; filename="./inpainted_ecco2_fields.nc", kw...)

    # Fields initialized from ECCO2
    grid = field.grid
    name = ecco2_metadata.name

    mask = ecco2_center_mask(architecture(grid))
    
    f = inpainted_ecco2_field(name; filename, mask,
                              architecture = architecture(grid),
                              kw...)

    f_grid = Field(ecco2_location[name], grid)   
    three_dimensional_interpolate!(f_grid, f)
    set!(field, f_grid)

    return field
end

end # module



