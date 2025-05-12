using NCDatasets
using JLD2
using ClimaOcean.InitialConditions: interpolate!

import Oceananigans.Fields: set!, Field

restrict(::Nothing, interfaces, N) = interfaces, N

# TODO support stretched native grids
function restrict(bbox_interfaces, interfaces, N)
    Δ = interfaces[2] - interfaces[1]
    rΔ = bbox_interfaces[2] - bbox_interfaces[1]
    ϵ = rΔ / Δ
    rN = Integer(ϵ * N)
    return bbox_interfaces, rN
end

function native_grid(metadata::Metadata, arch=CPU(); halo = (3, 3, 3))
    Nx, Ny, Nz = size(metadata)
    FT = eltype(metadata)

    longitude = longitude_interfaces(metadata)
    latitude = latitude_interfaces(metadata)
    z = z_interfaces(metadata)

    # Restrict with BoundingBox
    bbox = metadata.bounding_box
    if !isnothing(bbox)
        longitude, Nx = restrict(bbox.longitude, longitude, Nx)
        latitude, Ny = restrict(bbox.latitude, latitude, Ny)
        # TODO: restrict in z too
    end

    grid = LatitudeLongitudeGrid(arch, FT; halo, longitude, latitude, z,
                                 size = (Nx, Ny, Nz))

    return grid
end

"""
    Field(metadata::Metadatum;
          architecture = CPU(),
          inpainting = default_inpainting(metadata),
          mask = nothing,
          halo = (7, 7, 7),
          cache_inpainted_data = true)

Return a `Field` on `architecture` described by `metadata` with `halo` size.
If not `nothing`, the `inpainting` method is used to fill the cells
within the specified `mask`. `mask` is set to `compute_mask` for non-nothing
`inpainting`. Keyword argument `cache_inpainted_data` dictates whether the inpainted
data is cached to avoid recomputing it; default: `true`.
"""
function Field(metadata::Metadatum, arch=CPU();
               inpainting = default_inpainting(metadata),
               mask = nothing,
               halo = (3, 3, 3),
               cache_inpainted_data = true)

    grid = native_grid(metadata, arch; halo)
    LX, LY, LZ = location(metadata)
    field = Field{LX, LY, LZ}(grid)

    if !isnothing(inpainting)
        inpainted_path = inpainted_metadata_path(metadata)
        if isfile(inpainted_path)
            file = jldopen(inpainted_path, "r")
            maxiter = file["inpainting_maxiter"]

            # read data if generated with the same inpainting
            if maxiter == inpainting.maxiter
                data = file["data"]
                close(file)
                try
                    copyto!(parent(field), data)
                    return field
                catch
                    @warn "Could not load existing inpainted data at $inpainted_path.\n" *
                          "Re-inpainting and saving data..."
                    rm(inpainted_path, force=true)
                end
            end

            close(file)
        end
    end

    download_dataset(metadata)
    path = metadata_path(metadata)
    dsname = dataset_variable_name(metadata)

    # NetCDF shenanigans
    ds = Dataset(path)
    if is_three_dimensional(metadata)
        data = ds[dsname][:, :, :, 1]

        # Many ocean datasets use a "depth convention" for their vertical axis
        if reversed_vertical_axis(metadata.dataset)
            data = reverse(data, dims=3)
        end
    else
        data = ds[dsname][:, :, 1]
    end

    close(ds)

    set_metadata_field!(field, data, metadata)
    fill_halo_regions!(field)

    if !isnothing(inpainting)
        # Respect user-supplied mask, but otherwise build default mask for this dataset.
        if isnothing(mask)
            mask = compute_mask(metadata, field)
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

function set!(target_field::Field, metadata::Metadatum; kw...)
    grid = target_field.grid
    arch = child_architecture(grid)
    meta_field = Field(metadata, arch; kw...)
    interpolate!(target_field, meta_field)
    return target_field
end

# manglings
struct ShiftSouth end
struct AverageNorthSouth end

reversed_sign(dataset, val_name) = false

@inline mangle(i, j, data, ::Nothing) = @inbounds data[i, j]
@inline mangle(i, j, data, ::ShiftSouth) = @inbounds data[i, j-1]
@inline mangle(i, j, data, ::AverageNorthSouth) = @inbounds (data[i, j+1] + data[i, j]) / 2

@inline mangle(i, j, k, data, ::Nothing) = @inbounds data[i, j, k]
@inline mangle(i, j, k, data, ::ShiftSouth) = @inbounds data[i, j-1, k]
@inline mangle(i, j, k, data, ::AverageNorthSouth) = @inbounds (data[i, j+1, k] + data[i, j, k]) / 2

@inline maybe_reverse_sign(datum, reverse::Bool) = ifelse(reverse, - datum, datum)

function set_metadata_field!(field, data, metadatum)
    grid = field.grid
    arch = architecture(grid)

    Nx, Ny, Nz = size(metadatum)
    mangling = if size(data, 2) == Ny-1
        ShiftSouth()
    elseif size(data, 2) == Ny+1
        AverageNorthSouth()
    else
        nothing
    end

    temp_units = if metadatum.name == :temperature
        temperature_units(metadatum.dataset)
    else
        nothing
    end

    reverse_sign = reversed_sign(metadatum.dataset, Val(metadatum.name))

    if ndims(data) == 2
        _kernel = _set_2d_metadata_field!
        spec = :xy
    else
        _kernel = _set_3d_metadata_field!
        spec = :xyz
    end

    data = on_architecture(arch, data)
    Oceananigans.Utils.launch!(arch, grid, spec, _kernel, field, data, mangling, reverse_sign, temp_units)

    return nothing
end

@kernel function _set_2d_metadata_field!(field, data, mangling, reverse_sign, temp_units)
    i, j = @index(Global, NTuple)
    d = mangle(i, j, data, mangling)
    d = convert_temperature(d, temp_units)
    d = maybe_reverse_sign(d, reverse_sign)
    @inbounds field[i, j, 1] = d
end

@inline nan_convert_missing(FT, ::Missing) = convert(FT, NaN)
@inline nan_convert_missing(FT, d::Number) = convert(FT, d)

@kernel function _set_3d_metadata_field!(field, data, mangling, reverse_sign, temp_units)
    i, j, k = @index(Global, NTuple)
    FT = eltype(field)
    d = mangle(i, j, k, data, mangling)
    d = nan_convert_missing(FT, d)
    d = convert_temperature(d, temp_units)
    d = maybe_reverse_sign(d, reverse_sign)
    @inbounds field[i, j, k] = d
end

@inline convert_temperature(T, units) = T
@inline function convert_temperature(T::FT, ::Kelvin) where FT
    T₀ = convert(FT, 273.15)
    return T - T₀
end


#####
##### Masking data for inpainting
#####

"""
    compute_mask(metadata::Metadatum, dataset_field,
                 mask_value = default_mask_value(metadata),
                 minimum_value = -1f5,
                 maximum_value = 1f5)

A boolean field where `true` represents a missing value in the dataset_field.
"""
function compute_mask(metadata::Metadatum, dataset_field,
                      mask_value = default_mask_value(metadata.dataset),
                      minimum_value = -1f5,
                      maximum_value = 1f5)

    grid = dataset_field.grid
    arch = Oceananigans.Architectures.architecture(grid)
    LX, LY, LZ = location(dataset_field)
    mask = Field{LX, LY, LZ}(grid, Bool)

    # Set the mask with zeros where field is defined
    launch!(arch, grid, :xyz, _compute_mask!,
            mask, dataset_field, minimum_value, maximum_value, mask_value)

    return mask
end

@kernel function _compute_mask!(mask, field, min_value, max_value, mask_value)
    i, j, k = @index(Global, NTuple)
    @inbounds mask[i, j, k] = is_masked(field[i, j, k], min_value, max_value, mask_value)
end

@inline is_masked(a, min_value, max_value, mask_value) =
    isnan(a) | (a <= min_value) | (a >= max_value) | (a == mask_value)