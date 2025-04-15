module EN4

export EN4Metadatum, EN4_field, EN4_mask, EN4_immersed_grid, adjusted_EN4_tracers, initialize!
export EN4Monthly
export EN4FieldTimeSeries, EN4Restoring, LinearlyTaperedPolarMask

using ClimaOcean
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: inpaint_mask!, NearestNeighborInpainting,
                                download_progress, compute_native_date_range,
                                shift_longitude_to_0_360, Kelvin, Celsius

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

import ..z_faces, ..empty_field,
       ..variable_is_three_dimensional

download_EN4_cache::String = ""
function __init__()
    global download_EN4_cache = @get_scratch!("EN4")
end

include("EN4_metadata.jl")
include("EN4_mask.jl")

vertical_interfaces(metadata::Metadata{<:EN4Monthly}) =
    [
    -5500.0,
    -5200.5986,
    -4901.459,
    -4602.695,
    -4304.469,
    -4007.0154,
    -3710.6636,
    -3415.8828,
    -3123.3313,
    -2833.926,
    -2548.9243,
    -2270.0168,
    -1999.4135,
    -1739.8875,
    -1494.7288,
    -1267.542,
    -1061.846,
    -880.5102,
    -725.1709,
    -595.8501,
    -490.9494,
    -407.6244,
    -342.3651,
    -291.5638,
    -251.9186,
    -220.6418,
    -195.5097,
    -174.8189,
    -157.3019,
    -142.0345,
    -128.3522,
    -115.7825,
    -103.9913,
    -92.7437,
    -81.8753,
    -71.2711,
    -60.8509,
    -50.5586,
    -40.3553,
    -30.214,
    -20.1158,
    -10.0475,
    -0.0,
    ]

empty_EN4_field(variable_name::Symbol; kw...) = empty_field(Metadatum(variable_name, dataset=EN4Monthly()); kw...)

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
    EN4_field(metadata::EN4Metadatum;
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
function EN4_field(metadata::EN4Metadatum;
                   architecture = CPU(),
                   inpainting = default_inpainting(metadata),
                   mask = nothing,
                   horizontal_halo = (7, 7),
                   cache_inpainted_data = true)

    field = empty_field(metadata; architecture, horizontal_halo)
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

    # Convert data from Union{Missing, FT} to FT
    FT = eltype(field)
    data[ismissing.(data)] .= 1e10 # Artificially large number!
    data = if location(field)[2] == Face # ?
        new_data = zeros(FT, size(field))
        new_data[:, 1:end-1, :] .= data
        new_data
    else
        Array{FT}(data)
    end

    if metadata.name == :temperature && dataset_temperature_units(metadata) isa Kelvin
        data[data .!= FT(1e10)] .-= FT(273.15) # convert to Celsius
    end

    data = shift_longitude_to_0_360(data, metadata)

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
    var = string(metadata.name)
    return without_extension * "_" * var *"_inpainted.jld2"
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
