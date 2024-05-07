module ECCO

export ECCOMetadata, ecco_field, ecco_center_mask, adjusted_ecco_tracers, initialize!

using ClimaOcean.DataWrangling: inpaint_mask!
using ClimaOcean.InitialConditions: three_dimensional_regrid!, interpolate!

using Oceananigans
using Oceananigans.Architectures: architecture, child_architecture
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations: DistributedField, all_reduce, barrier!
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index
using NCDatasets
using Downloads: download
using Dates

using Adapt

include("ecco_metadata.jl")

# Vertical coordinate
const ECCO_z = [
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

empty_ecco_field(variable_name::Symbol; kw...) = empty_ecco_field(ECCOMetadata(variable_name); kw...)

function empty_ecco_field(metadata::ECCOMetadata;
                          architecture = CPU(), 
                          horizontal_halo = (3, 3))

    Nx, Ny, Nz, _ = size(metadata)

    variable_name = metadata.name
    location = field_location(metadata)
    
    location = ecco_location[variable_name]

    longitude = (0, 360)
    latitude = (-90, 90)
    TX, TY = (Periodic, Bounded)

    if variable_is_three_dimensional(metadata)
        z    = ECCO_z
        # add vertical halo for 3D fields
        halo = (horizontal_halo..., 3)
        LZ   = Center
        TZ   = Bounded
        N    = (Nx, Ny, Nz)
    else
        z    = nothing
        halo = horizontal_halo
        LZ   = Nothing
        TZ   = Flat
        N    = (Nx, Ny)
    end

    # Flat in z if the variable is two-dimensional
    grid = LatitudeLongitudeGrid(architecture; halo, size = N, topology = (TX, TY, TZ),
                                 longitude, latitude, z)

    return Field{location...}(grid)
end

"""
    ecco_field(variable_name;
                architecture = CPU(),
                horizontal_halo = (1, 1),
                user_data = nothing,
                url = ecco_urls[variable_name],
                filename = ecco_file_names[variable_name],
                short_name = ecco_short_names[variable_name])

Retrieve the ecco field corresponding to `variable_name`. 
The data is either:
(1) retrieved from `filename`,
(2) dowloaded from `url` if `filename` does not exists,
(3) filled from `user_data` if `user_data` is provided.
"""
function ecco_field(metadata::ECCOMetadata;
                     architecture = CPU(),
                     horizontal_halo = (3, 3),
                     filename = file_name(metadata))

    shortname = short_name(metadata)
    
    download_dataset!(metadata)

    ds = Dataset(filename)
    if variable_is_three_dimensional(metadata)
        data = ds[shortname][:, :, :, 1]
        # The surface layer in three-dimensional ECCO fields is at `k = 1`
        data = reverse(data, dims = 3)
    else
        data = ds[shortname][:, :, 1]
    end        
    close(ds)

    field = empty_ecco_field(metadata; architecture, horizontal_halo)
    
    FT    = eltype(field)
    data[ismissing.(data)] .= 1e10 # Artificially large number!
    data  = if location(field)[2] == Face
        new_data = zeros(FT, size(field))
        new_data[:, 1:end-1, :] .= data
        new_data    
    else
        convert.(FT, data)
    end

    set!(field, data)
    fill_halo_regions!(field)

    return field
end

@kernel function _set_ecco_mask!(mask, Tᵢ, minimum_value, maximum_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = (Tᵢ[i, j, k] < minimum_value) | (Tᵢ[i, j, k] > maximum_value) 
end

"""
    ecco_center_mask(architecture = CPU(); minimum_value = Float32(-1e5))

A boolean field where `false` represents a missing value in the ECCO dataset.
"""
function ecco_center_mask(architecture = CPU(); 
                          minimum_value = Float32(-1e5),
                          maximum_value = Float32(1e5),
                          metadata = ECCOMetadata(:temperature),
                          filename = file_name(metadata))

    field = ecco_field(metadata; architecture, filename)
    mask  = CenterField(field.grid, Bool)

    # Set the mask with ones where T is defined
    launch!(architecture, field.grid, :xyz, _set_ecco_mask!, mask, field, minimum_value, maximum_value)

    return mask
end

"""
    inpainted_ecco_field(variable_name; 
                          architecture = CPU(),
                          filename = "./inpainted_ecco_fields.nc",
                          mask = ecco_center_mask(architecture))
    
Retrieve the ECCO field corresponding to `variable_name` inpainted to fill all the
missing values in the original dataset.

Arguments:
==========

- `variable_name`: the variable name corresponding to the Dataset.

Keyword Arguments:
==================

- `architecture`: either `CPU()` or `GPU()`.

- `filename`: the path where to retrieve the data from. If the file does not exist,
              the data will be downloaded from the ECCO dataset.

- `mask`: the mask used to inpaint the field (see `inpaint_mask!`).

- `maxiter`: the maximum number of iterations to inpaint the field (see `inpaint_mask!`).

"""
function inpainted_ecco_field(metadata::ECCOMetadata; 
                              architecture = CPU(),
                              filename = file_name(metadata),
                              mask = ecco_center_mask(architecture),
                              maxiter = 10,
                              kw...)
    
    f = ecco_field(metadata; architecture, filename, kw...)
        
    # Make sure all values are extended properly
    @info "In-painting ecco $(metadata.name)"
    inpaint_mask!(f, mask; maxiter)

    fill_halo_regions!(f)

    return f
end

inpainted_ecco_field(variable_name::Symbol; kw...) = inpainted_ecco_field(ECCOMetadata(variable_name); kw...)
    
function set!(field::DistributedField, ecco_metadata::ECCOMetadata; kw...)
    # Fields initialized from ECCO
    grid = field.grid
    arch = architecture(grid)
    child_arch = child_architecture(arch)

    f_ecco = if arch.local_rank == 0 # Make sure we read/write the file using only one core
        mask = ecco_center_mask(child_arch)
        
        inpainted_ecco_field(ecco_metadata; mask,
                                  architecture = child_arch,
                                  kw...)
    else
        empty_ecco_field(ecco_metadata; architecture = child_arch)
    end

    barrier!(arch)

    # Distribute ecco field to all workers
    parent(f_ecco) .= all_reduce(+, parent(f_ecco), arch)

    f_grid = Field(field_location(ecco_metadata), grid)   
    interpolate!(f_grid, f_ecco)
    set!(field, f_grid)
    
    return field
end

function set!(field::Field, ecco_metadata::ECCOMetadata; kw...)

    # Fields initialized from ECCO
    grid = field.grid
    mask = ecco_center_mask(architecture(grid))
    
    f = inpainted_ecco_field(ecco_metadata; mask,
                              architecture = architecture(grid),
                              kw...)

    f_grid = Field(field_location(ecco_metadata), grid)   
    interpolate!(f_grid, f)
    set!(field, f_grid)

    return field
end

include("ecco_restoring.jl")

end # Module 