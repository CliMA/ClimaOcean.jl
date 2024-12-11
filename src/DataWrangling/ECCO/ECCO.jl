module ECCO

export ECCOMetadata, ECCO_field, ECCO_mask, ECCO_immersed_grid, adjusted_ECCO_tracers, initialize!
export ECCO2Monthly, ECCO4Monthly, ECCO2Daily
export ECCORestoring, LinearlyTaperedPolarMask

using ClimaOcean
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: inpaint_mask!, NearestNeighborInpainting, download_progress
using ClimaOcean.InitialConditions: three_dimensional_regrid!, interpolate!

using Oceananigans
using Oceananigans: location
using Oceananigans.Architectures: architecture, child_architecture
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations: DistributedField, all_reduce, barrier!
using Oceananigans.Utils

using KernelAbstractions: @kernel, @index
using NCDatasets
using JLD2
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
                          horizontal_halo = (7, 7))

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
    ECCO_field(metadata::ECCOMetadata;
               architecture = CPU(),
               inpainting = nothing,
               mask = nothing,
               horizontal_halo = (7, 7),
               cache_inpainted_data = false)

Return a `Field` on `architecture` described by `ECCOMetadata` with
`horizontal_halo` size.
If not `nothing`, the `inpainting` method is used to fill the cells
within the specified `mask`. `mask` is set to `ECCO_mask` for non-nothing
`inpainting`.
"""
function ECCO_field(metadata::ECCOMetadata;
                    architecture = CPU(),
                    inpainting = NearestNeighborInpainting(Inf),
                    mask = nothing,
                    horizontal_halo = (7, 7),
                    cache_inpainted_data = true)

    field = empty_ECCO_field(metadata; architecture, horizontal_halo)
    inpainted_path = inpainted_metadata_path(metadata)

    if !isnothing(inpainting) && isfile(inpainted_path)
        file = jldopen(inpainted_path, "r")
        maxiter = file["inpainting_maxiter"]

        # read data if generated with the same inpainting
        if maxiter == inpainting.maxiter
            data = file["data"]
            close(file)
            copyto!(parent(field), data)
            return field
        end

        close(file)
    end

    download_dataset(metadata)
    path = metadata_path(metadata)
    ds = Dataset(path)
    shortname = short_name(metadata)

    if variable_is_three_dimensional(metadata)
        data = ds[shortname][:, :, :, 1]
        data = reverse(data, dims=3)
    else
        data = ds[shortname][:, :, 1]
    end        

    close(ds)
    
    # Convert data from Union(FT, missing} to FT
    FT = eltype(field)
    data[ismissing.(data)] .= 1e10 # Artificially large number!
    data = if location(field)[2] == Face # ?
        new_data = zeros(FT, size(field))
        new_data[:, 1:end-1, :] .= data
        new_data    
    else
        data = Array{FT}(data)
    end
    
    # ECCO4 data is on a -180, 180 longitude grid as opposed to ECCO2 data that
    # is on a 0, 360 longitude grid. To make the data consistent, we shift ECCO4
    # data by 180 degrees in longitude
    if metadata.version isa ECCO4Monthly 
        Nx = size(data, 1)
        if variable_is_three_dimensional(metadata)
            shift = (Nx ÷ 2, 0, 0)
        else
            shift = (Nx ÷ 2, 0)
        end
        data = circshift(data, shift)
    end

    set!(field, data)
    fill_halo_regions!(field)

    if !isnothing(inpainting)
        # Respect user-supplied mask, but otherwise build default ECCO mask.
        if isnothing(mask)
            mask = ECCO_mask(metadata, architecture; data_field=field)
        end

        # Make sure all values are extended properly
        name = string(metadata.name)
        date = string(metadata.dates)
        version = summary(metadata.version)
        @info string("Inpainting ", version, " ", name, " data from ", date, "...")
        start_time = time_ns()
        
        inpaint_mask!(field, mask; inpainting)
        fill_halo_regions!(field)

        elapsed = 1e-9 * (time_ns() - start_time)
        @info string(" ... (", prettytime(elapsed), ")")
    
        # We cache the inpainted data to avoid recomputing it
        @root if cache_inpainted_data
            file = jldopen(inpainted_path, "w+")
            file["data"] = on_architecture(CPU(), parent(field))
            file["inpainting_maxiter"] = inpainting.maxiter
            close(file)
        end
    end

    return field
end

# Fallback
ECCO_field(var_name::Symbol; kw...) = ECCO_field(ECCOMetadata(var_name); kw...)

function inpainted_metadata_filename(metadata::ECCOMetadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::ECCOMetadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

function set!(field::Field, ECCO_metadata::ECCOMetadata; kw...)

    # Fields initialized from ECCO
    grid = field.grid
    arch = child_architecture(grid)
    mask = ECCO_mask(ECCO_metadata, arch)

    f = ECCO_field(ECCO_metadata; mask,
                   architecture = arch,
                   kw...)

    interpolate!(field, f)

    return field
end

end # Module 

