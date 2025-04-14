module EN4

export EN4Metadatum, EN4_field, EN4_mask, EN4_immersed_grid, adjusted_EN4_tracers, initialize!
export EN4Monthly
export EN4FieldTimeSeries, EN4Restoring, LinearlyTaperedPolarMask

using ClimaOcean
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: inpaint_mask!, NearestNeighborInpainting, download_progress, compute_native_date_range
using ClimaOcean.InitialConditions: three_dimensional_regrid!, interpolate!

using Oceananigans
using Oceananigans: location
using Oceananigans.Architectures: architecture, child_architecture
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: DistributedField, all_reduce, barrier!
using Oceananigans.Utils

using KernelAbstractions: @kernel, @index
using NCDatasets
using JLD2
using Downloads: download
using Dates
using Adapt
using Scratch

download_EN4_cache::String = ""
function __init__()
    global download_EN4_cache = @get_scratch!("EN4")
end

include("EN4_metadata.jl")
include("EN4_mask.jl")

# Vertical coordinate
const EN4_z = [
    -5350.272,
    -5050.9897,
    -4752.021,
    -4453.502,
    -4155.628,
    -3858.6763,
    -3563.0408,
    -3269.278,
    -2978.166,
    -2690.7803,
    -2408.5835,
    -2133.517,
    -1868.0707,
    -1615.2905,
    -1378.661,
    -1161.8059,
    -967.99585,
    -799.5496,
    -657.32294,
    -540.5022,
    -446.80093,
    -372.96545,
    -315.37408,
    -270.53412,
    -235.38617,
    -207.42545,
    -184.69746,
    -165.72845,
    -149.43373,
    -135.02855,
    -121.9519,
    -109.806175,
    -98.31118,
    -87.27029,
    -76.54591,
    -66.041985,
    -55.691494,
    -45.44776,
    -35.27829,
    -25.16046,
    -15.07854,
    -5.0215898,
]

empty_EN4_field(variable_name::Symbol; kw...) = empty_EN4_field(Metadatum(variable_name, dataset=EN4Monthly()); kw...)

function empty_EN4_field(metadata::EN4Metadata;
                          architecture = CPU(), 
                          horizontal_halo = (7, 7))

    Nx, Ny, Nz, _ = size(metadata)
    loc = location(metadata)
    longitude = (0, 360)
    latitude = (-83, 89)
    TX, TY = (Periodic, Bounded)

    if variable_is_three_dimensional(metadata)
        TZ = Bounded
        LZ = Center
        z = EN4_z
        halo = (horizontal_halo..., 3)
        sz = (Nx, Ny, Nz)
    else # the variable is two-dimensional
        TZ = Flat
        LZ = Nothing
        z = nothing
        halo = horizontal_halo
        sz = (Nx, Ny)
    end

    grid = LatitudeLongitudeGrid(architecture, Float32; halo, longitude, latitude, z,
                                 size = sz,
                                 topology = (TX, TY, TZ))

    return Field{loc...}(grid)
end

# Only temperature and salinity need a thorough inpainting because of stability,
# other variables can do with only a couple of passes. Sea ice variables 
# cannot be inpainted because zeros in the data are physical, not missing values.
function default_inpainting(metadata::EN4Metadata)
    if metadata.name in [:temperature, :salinity]
        return NearestNeighborInpainting(Inf)
    elseif metadata.name in [:sea_ice_fraction, :sea_ice_thickness]
        return nothing
    else
        return NearestNeighborInpainting(5)
    end
end

"""
    EN4_field(metadata::EN4Metadata;
               architecture = CPU(),
               inpainting = nothing,
               mask = nothing,
               horizontal_halo = (7, 7),
               cache_inpainted_data = false)

Return a `Field` on `architecture` described by `EN4Metadata` with
`horizontal_halo` size.
If not `nothing`, the `inpainting` method is used to fill the cells
within the specified `mask`. `mask` is set to `EN4_mask` for non-nothing
`inpainting`.
"""
function EN4_field(metadata::EN4Metadata;
                    architecture = CPU(),
                    inpainting = default_inpainting(metadata),
                    mask = nothing,
                    horizontal_halo = (7, 7),
                    cache_inpainted_data = true)
                    
    field = empty_EN4_field(metadata; architecture, horizontal_halo)
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
    
    set!(field, data)
    fill_halo_regions!(field)

    if !isnothing(inpainting)
        # Respect user-supplied mask, but otherwise build default EN4 mask.
        if isnothing(mask)
            mask = EN4_mask(metadata, architecture; data_field=field)
        end

        # Make sure all values are extended properly
        name = string(metadata.name)
        date = string(metadata.dates)
        dataset = summary(metadata.dataset)
        @info string("Inpainting ", dataset, " ", name, " data from ", date, "...")
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
EN4_field(var_name::Symbol; kw...) = EN4_field(EN4Metadata(var_name); kw...)

function inpainted_metadata_filename(metadata::EN4Metadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::EN4Metadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

function set!(field::Field, EN4_metadata::EN4Metadatum; kw...)

    # Fields initialized from EN4
    grid = field.grid
    arch = child_architecture(grid)
    mask = EN4_mask(EN4_metadata, arch)

    f = EN4_field(EN4_metadata; mask,
                   architecture = arch,
                   kw...)

    interpolate!(field, f)

    return field
end

include("EN4_restoring.jl")

end # Module 

