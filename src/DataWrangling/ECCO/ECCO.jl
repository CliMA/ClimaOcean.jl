module ECCO

export ECCOMetadata, ecco_field, ECCO_missings_field, adjusted_ecco_tracers, initialize!
export ECCO2Monthly, ECCO4Monthly, ECCO2Daily
export ECCO_restoring_forcing

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
include("ECCO_missings_field.jl")

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
    location = location(metadata)
    
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
                filename = ecco_metadata_filenames[variable_name],
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
                    filename = metadata_filename(metadata))

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
    
    # ECCO4 data is on a -180, 180 longitude grid as opposed to ECCO2 data that
    # is on a 0, 360 longitude grid. To make the data consistent, we shift ECCO4
    # data by 180 degrees in longitude
    if metadata.version isa ECCO4Monthly 
        Nx = size(data, 1)
        data = circshift(data, (Nx ÷ 2, 0, 0))
    end

    set!(field, data)
    fill_halo_regions!(field)

    return field
end

# Fallback
ecco_field(var_name::Symbol; kw...) = ecco_field(ECCOMetadata(var_name); kw...)

"""
    inpainted_ecco_field(variable_name; 
                          architecture = CPU(),
                          filename = "./inpainted_ecco_fields.nc",
                          mask = ECCO_missings_field(architecture))
    
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
                              filename = metadata_filename(metadata),
                              mask = ECCO_missings_field(metadata, architecture),
                              maxiter = Inf,
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

    ecco_field = if arch.local_rank == 0 # Make sure we read/write the file using only one core
        mask = ECCO_missings_field(ecco_metadata, child_arch)
        inpainted_ecco_field(ecco_metadata; mask, architecture=child_arch, kw...)
    else
        empty_ecco_field(ecco_metadata; architecture = child_arch)
    end

    barrier!(arch)

    # Distribute ecco field to all workers
    parent(ecco_field) .= all_reduce(+, parent(ecco_field), arch)

    ecco_loc = location(ecco_metadata)
    if 
    grid_field = Field(location(ecco_metadata), grid)   
    interpolate!(grid_field, ecco_field)
    set!(field, grid_field)
    
    return field
end

function set!(field::Field, metadata::ECCOMetadata; kw...)

    # Fields initialized from ECCO
    grid = field.grid
    mask = ECCO_missings_field(metadata, architecture(grid))
    ecco_field = inpainted_ecco_field(metadata; mask, architecture=architecture(grid), kw...)

    if location(field) === location(metadata)
        interpolate!(field, ecco_field)
    else # how does this help?
        field_loc = Field(location(metadata), grid)   
        interpolate!(field_loc, ecco_field)
        set!(field, field_loc)
    end

    return field
end

include("ecco_restoring.jl")

end # Module 
