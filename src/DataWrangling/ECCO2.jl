module ECCO2

export ECCO2Metadata, ecco2_field, ecco2_center_mask, adjusted_ecco_tracers, initialize!

using ClimaOcean.DataWrangling: fill_missing_values!
using ClimaOcean.InitialConditions: three_dimensional_regrid!

using Oceananigans
using Oceananigans: architecture
using Oceananigans.BoundaryConditions
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index
using NCDatasets

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

ecco2_tracer_fields = Dict(
    :ecco2_temperature => :temperature,
    :ecco2_salinity => :salinity,
    :ecco_2_effective_ice_thickness => :effective_ice_thickness
)

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

ecco2_location = Dict(
    :temperature   => (Center, Center, Center),
    :salinity      => (Center, Center, Center),
    :effective_ice_thickness => (Center, Center, Nothing)
)

ecco2_depth_names = Dict(
    :temperature   => "DEPTH_T",
    :salinity      => "DEPTH_T",
)

shortnames = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_area_fraction => "SIarea",
    :u_velocity            => "UVEL",
    :v_velocity            => "VVEL",
)


function empty_ecco2_field(data::ECCO2Metadata; architecture = CPU(), 
                                            horizontal_halo = (1, 1))

    variable_name = data.name

    location = ecco2_location[variable_name]

    grid = LatitudeLongitudeGrid(architecture,
                                 size = (ECCO2_Nx, ECCO2_Ny, ECCO2_Nz),
                                 longitude = (0, 360),
                                 latitude = (-90, 90),
                                 z = ECCO2_z,
                                 halo = (1, 1, 1),
                                 topology = (Periodic, Bounded, Bounded))

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

Retrieve the ecco2 field corresponding to `variable_name`. The data is either:
(1) retrieved from `filename` 
(2) dowloaded from `url` if `filename` does not exists
(3) filled from `user_data` if `user_data` is provided
"""
function ecco2_field(variable_name;
                     architecture = CPU(),
                     horizontal_halo = (1, 1),
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
    @inbounds mask[i, j, k] = ifelse(Tᵢ[i, j, k] < minimum_value, 0, 1)
end

"""
    ecco2_center_mask(architecture = CPU(); minimum_value = Float32(-1e5))

An integer field where 0 represents a missing value in the ECCO2 :temperature
dataset and 1 represents a valid value
"""
function ecco2_center_mask(architecture = CPU(); minimum_value = Float32(-1e5))
    Tᵢ   = ecco2_field(:temperature; architecture)
    mask = CenterField(Tᵢ.grid)

    # Set the mask with ones where T is defined
    launch!(architecture, Tᵢ.grid, :xyz, _set_ecco2_mask!, mask, Tᵢ, minimum_value)

    return mask
end

"""
    adjusted_ecco_field(variable_name; 
                        architecture = CPU(),
                        filename = "./data/initial_ecco_tracers.nc",
                        mask = ecco2_center_mask(architecture))
    
Retrieve the ECCO2 field corresponding to `variable_name` adjusted to fill all the
missing values in the original dataset

Arguments:
==========

- `variable_name`: the variable name corresponding to the Dataset

Keyword Arguments:
==================

- `architecture`: either `CPU()` or `GPU()`

- `filename`: the path where to retrieve the data from. If the file does not exist,
              the data will be retrived from the ECCO2 dataset, it will be adjusted and
              saved down in `filename`

- `mask`: the mask used to extend the field (see `adjust_tracer!`)
"""
function adjusted_ecco2_field(variable_name; 
                              architecture = CPU(),
                              filename = "./data/adjusted_ecco_tracers.nc",
                              mask = ecco2_center_mask(architecture))
    
    if !isfile(filename)
        f = ecco2_field(variable_name; architecture)
        
        # Make sure all values are extended properly
        @info "in-painting ecco field $variable_name and saving it in $filename"
        fill_missing_values!(f; mask)

        ds = Dataset(filename, "c")
        defVar(ds, string(variable_name), Array(interior(f)), ("lat", "lon", "z"))

        close(ds)
    else
        ds = Dataset(filename, "a")

        if haskey(ds, string(variable_name))
            data = ds[variable_name][:, :, :]
            f = ecco2_field(variable_name; architecture, user_data = data)
        else
            f = ecco2_field(variable_name; architecture)
            # Make sure all values are extended properly
            @info "in-painting ecco field $variable_name and saving it in $filename"
            fill_missing_values!(f; mask)

            defVar(ds, string(variable_name), Array(interior(f)), ("lat", "lon", "z"))
        end

        close(ds)
    end

    return f
end

function set!(field::Field, ecco2::ECCO2Metadata; filename = "./data/adjusted_ecco_tracers.nc")
    # Fields initialized from ECCO2
    grid = field.grid

    mask = ecco2_center_mask(architecture(grid))
    
    f = adjusted_ecco2_field(ecco2.name; 
                             architecture = architecture(grid),
                             filename,
                             mask)

    f_grid = Field(ecco2_location[ecco2.name], grid)   

    three_dimensional_regrid!(f_grid, f)

    set!(field, f_grid)

    return field
end

end # module

