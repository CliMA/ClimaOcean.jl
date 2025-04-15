using NCDatasets
using JLD2

function inpainted_metadata_path end
function default_inpainting end

"""
    dataset_field(metadata::Metadatum;
                  architecture = CPU(),
                  inpainting = nothing,
                  mask = nothing,
                  horizontal_halo = (7, 7),
                  cache_inpainted_data = false)

Return a `Field` on `architecture` described by `metadata` with
`horizontal_halo` size.
If not `nothing`, the `inpainting` method is used to fill the cells
within the specified `mask`. `mask` is set to `Dataset_mask` for non-nothing
`inpainting`.
"""
function dataset_field(metadata::Metadatum;
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
