using ClimaOcean.InitialConditions: interpolate!

using NCDatasets
using JLD2

import Oceananigans.Fields: set!, Field

function inpainted_metadata_path end

""" Amount to shift the data in longitude (in degrees East). """
longitude_shift(metadata) = 0

function shift_longitude_to_0_360(data, metadata)
    Nx = size(data, 1)

    degrees = longitude_shift(metadata)
    Nshift = Integer(degrees / 360 * Nx)

    if variable_is_three_dimensional(metadata)
        shift = (Nshift, 0, 0)
    else
        shift = (Nshift, 0)
    end
    data = circshift(data, shift)

    return data
end

# Only temperature and salinity need a thorough inpainting because of stability,
# other variables can do with only a couple of passes. Sea ice variables
# cannot be inpainted because zeros in the data are physical, not missing values.
function default_inpainting(metadata)
    if metadata.name in [:temperature, :salinity]
        return NearestNeighborInpainting(Inf)
    elseif metadata.name in [:sea_ice_thickness, :sea_ice_concentration]
        return nothing
    else
        return NearestNeighborInpainting(5)
    end
end

# Will be extended by the particular datase to
# make sure the data is in the right format and follows the
# ClimaOcean conventions.
enforce_data_conventions(data, dataset, val_name) = data

"""
    Field(metadata::Metadatum;
          architecture = CPU(),
          inpainting = default_inpainting(metadata),
          mask = nothing,
          horizontal_halo = (7, 7),
          cache_inpainted_data = true)

Return a `Field` on `architecture` described by `metadata` with `horizontal_halo` size.
If not `nothing`, the `inpainting` method is used to fill the cells
within the specified `mask`. `mask` is set to `dataset_mask` for non-nothing
`inpainting`. Keyword argument `cache_inpainted_data` dictates whether the inpainted
data is cached to avoid recomputing it; default: `true`.
"""
function Field(metadata::Metadatum;
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

    data = enforce_data_conventions(data, metadata.dataset, Val(metadata.name))
    data = shift_longitude_to_0_360(data, metadata)

    set!(field, data)
    fill_halo_regions!(field)

    if !isnothing(inpainting)
        # Respect user-supplied mask, but otherwise build default mask for this dataset.
        if isnothing(mask)
            mask = dataset_mask(metadata, architecture; data_field=field)
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

function set!(field::Field, metadata::Metadatum; kw...)

    # Fields initialized from metadata.dataset
    grid = field.grid
    arch = child_architecture(grid)
    mask = dataset_mask(metadata, arch)

    f = Field(metadata; mask,
              architecture = arch,
              kw...)

    interpolate!(field, f)

    return field
end
