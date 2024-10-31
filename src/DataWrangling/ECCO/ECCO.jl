module ECCO

export ECCOMetadata, ECCO_field, ECCO_mask, adjusted_ECCO_tracers, initialize!
export ECCO2Monthly, ECCO4Monthly, ECCO2Daily
export ECCO_restoring_forcing

using ClimaOcean
using ClimaOcean.DataWrangling
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
using Scratch

download_ECCO_cache::String = ""
function __init__()
    global download_ECCO_cache = @get_scratch!("ECCO")
end

include("ECCO_metadata.jl")
include("ECCO_mask.jl")
include("ECCO_restoring.jl")

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

empty_ECCO_field(variable_name::Symbol; kw...) = empty_ECCO_field(ECCOMetadata(variable_name); kw...)

function empty_ECCO_field(metadata::ECCOMetadata;
                          architecture = CPU(), 
                          horizontal_halo = (3, 3))

    Nx, Ny, Nz, _ = size(metadata)
    loc = location(metadata)
    longitude = (0, 360)
    latitude = (-90, 90)
    TX, TY = (Periodic, Bounded)

    if variable_is_three_dimensional(metadata)
        TZ = Bounded
        LZ = Center
        z = ECCO_z
        halo = (horizontal_halo..., 3)
        sz = (Nx, Ny, Nz)
    else # the variable is two-dimensional
        TZ = Flat
        LZ = Nothing
        z = nothing
        halo = horizontal_halo
        sz = (Nx, Ny)
    end

    grid = LatitudeLongitudeGrid(architecture; halo, longitude, latitude, z,
                                 size = sz,
                                 topology = (TX, TY, TZ))

    return Field{loc...}(grid)
end

"""
    ECCO_field(variable_name;
               architecture = CPU(),
               horizontal_halo = (1, 1),
               user_data = nothing,
               url = ecco_urls[variable_name],
               short_name = ecco_short_names[variable_name])

Retrieve the ECCO field corresponding to `variable_name`. 
The data is either:
(1) retrieved from `filename`,
(2) dowloaded from `url` if `filename` does not exists,
(3) filled from `user_data` if `user_data` is provided.
"""
function ECCO_field(metadata::ECCOMetadata;
                    architecture = CPU(),
                    horizontal_halo = (3, 3))

    download_dataset(metadata)
    path = metadata_path(metadata)
    ds = Dataset(path)

    shortname = short_name(metadata)

    if variable_is_three_dimensional(metadata)
        data = ds[shortname][:, :, :, 1]

        # The surface layer in three-dimensional ECCO fields is at `k = 1`
        data = reverse(data, dims=3)
    else
        data = ds[shortname][:, :, 1]
    end        

    close(ds)

    field = empty_ECCO_field(metadata; architecture, horizontal_halo)
    
    FT = eltype(field)
    data[ismissing.(data)] .= 1e10 # Artificially large number!
    data = if location(field)[2] == Face
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
ECCO_field(var_name::Symbol; kw...) = ECCO_field(ECCOMetadata(var_name); kw...)

"""
    inpainted_ECCO_field(variable_name; 
                         architecture = CPU(),
                         mask = ECCO_mask(architecture),
                         maxiter = Inf)
    
Retrieve the ECCO field corresponding to `variable_name` inpainted to fill all the
missing values in the original dataset.

Arguments:
==========

- `variable_name`: the variable name corresponding to the Dataset.

Keyword Arguments:
==================

- `architecture`: either `CPU()` or `GPU()`.
- `mask`: the mask used to inpaint the field (see `inpaint_mask!`).
- `maxiter`: the maximum number of iterations to inpaint the field (see `inpaint_mask!`).

"""
function inpainted_ECCO_field(metadata::ECCOMetadata; 
                              architecture = CPU(),
                              mask = ECCO_mask(metadata, architecture),
                              maxiter = Inf,
                              kw...)
    
    f = ECCO_field(metadata; architecture, kw...)

    # Make sure all values are extended properly
    @info "In-painting ECCO $(metadata.name)"
    inpaint_mask!(f, mask; maxiter)

    fill_halo_regions!(f)
    
    return f
end

inpainted_ECCO_field(variable_name::Symbol; kw...) = inpainted_ECCO_field(ECCOMetadata(variable_name); kw...)
    
function set!(field::Field, ecco_metadata::ECCOMetadata; kw...)

    # Fields initialized from ECCO
    grid = field.grid
    arch = child_architecture(grid)
    mask = ECCO_mask(ecco_metadata, arch)

    f = inpainted_ECCO_field(ecco_metadata; mask,
                             architecture = arch,
                             kw...)

    interpolate!(field, f)

    return field
end

end # Module 
