using ClimaOcean.InitialConditions: interpolate!

using NCDatasets
using JLD2

import Oceananigans.Fields: set!, Field

function inpainted_metadata_path end
mangle_data(metadata, data) = data

"""
    reversed_vertical_axis(metadata)

Return `true` if the vertical axis is reversed, `false` otherwise.
"""
reversed_vertical_axis(metadata) = false

""" Amount to shift the data in longitude (in degrees East). """
longitude_shift(metadata) = 0

function shift_longitude_to_0_360(data, metadata)
    Nx = size(data, 1)

    degrees = longitude_shift(metadata)
    Nshift = Integer(degrees / 360 * Nx)

    if is_three_dimensional(metadata)
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
               halo = (7, 7, 7),
               cache_inpainted_data = true)

    if !isnothing(inpainting)
        inpainted_path = inpainted_metadata_path(metadata)
        if isfile(inpainted_path)
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
    end

    download_dataset(metadata)
    path = metadata_path(metadata)
    dsname = dataset_variable_name(metadata)

    # NetCDF shenanigans
    ds = Dataset(path)
    if is_three_dimensional(metadata)
        data = ds[dsname][:, :, :, 1]

        # Many ocean datasets use a "depth convention" for their vertical axis
        if reversed_vertical_axis(metadata)
            data = reverse(data, dims=3)
        end
    else
        data = ds[dsname][:, :, 1]
    end

    close(ds)

    # Convert data from Union{Missing, FT} to FT
    data[ismissing.(data)] .= NaN
    data = shift_longitude_to_0_360(data, metadata)
    field = empty_field(metadata, arch; halo)
    FT = eltype(metadata)

    #=
    # Mangle the data, sadly.
    Nx, Ny, Nz = size(metadatum)
    data = if size(data, 2) == Ny-1 # seems to happen for Face fields
        # new_data = zeros(FT, Nx, Ny, Nz)
        new_data[:, 1:end-1, :] .= data
        new_data
    elseif size(data, 2) == Ny+1 # this happens for Copernicus...
        new_data = zeros(FT, Nx, Ny, Nz)
        new_data .= (data[:, 1:end-1, :] .+ data[:, 2:end, :]) ./ 2
        new_data
    else
        Array{FT}(data)
    end

    FT_NaN = convert(FT, NaN)
    if metadata.name == :temperature && dataset_temperature_units(metadata) isa Kelvin
        data[data .!= FT_NaN] .-= FT(273.15) # convert to Celsius
    end

    set!(field, data)
    =#

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

function set!(field::Field, metadata::Metadatum; kw...)

    # Fields initialized from metadata.dataset
    grid = field.grid
    arch = child_architecture(grid)
    mask = compute_mask(metadata, arch)

    f = Field(metadata, arch; mask, kw...)
    interpolate!(field, f)

    return field
end

# manglings
struct ShiftSouth end
struct AverageNorthSouth end

mangle(i, j, data, ::Nothing) = @inbounds data[i, j]
mangle(i, j, data, ::ShiftSouth) = @inbounds data[i, j-1]
mangle(i, j, data, ::AverageNorthSouth) = @inbounds (data[i, j] + data[i, j-1]) / 2

mangle(i, j, k, data, ::Nothing) = @inbounds data[i, j, k]
mangle(i, j, k, data, ::ShiftSouth) = @inbounds data[i, j-1, k]
mangle(i, j, k, data, ::AverageNorthSouth) = @inbounds (data[i, j-1, k] + data[i, j, k]) / 2

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

    temperature_units = if metadatum.name == :temperature
        dataset_temperature_units(metadata)
    else
        nothing
    end

    if ndims(data) == 2
        _kernel = _set_2d_metadata_field!
        spec = :xy
    else
        _kernel = _set_3d_metadata_field!
        spec = :xyz
    end

    Oceananigans.Utils.launch!(arch, grid, spec, _kernel, field, data, mangling, temperature_units)

    return nothing
end

@kernel function _set_2d_metadata_field!(field, data, mangling, temperature_units)
    i, j = @index(Global, NTuple)
    d = mangle(i, j, data, mangling)
    d = convert_temperature(d, temperature_units)
    @inbounds field[i, j, 1] = d
end

@kernel function _set_3d_metadata_field!(field, data, mangling, temperature_units)
    i, j, k = @index(Global, NTuple)
    d = mangle(i, j, k, data, mangling)
    d = convert_temperature(d, temperature_units)
    @inbounds field[i, j, k] = d
end

@inline convert_temperature(f, units) = f

@inline function convert_temperature(f::FT, ::Kelvin) where FT
    FTNaN = convert(FT, NaN)
    temperature_shift = convert(FT, 273.15)
    f = ifelse(f != FTNaN, f - temperature_shift, f)
    return f
end

#=
    # Mangle the data, sadly.
    Nx, Ny, Nz = size(metadatum)
    data = if size(data, 2) == Ny-1 # seems to happen for Face fields
        # new_data = zeros(FT, Nx, Ny, Nz)
        new_data[:, 1:end-1, :] .= data
        new_data
    elseif size(data, 2) == Ny+1 # this happens for Copernicus...
        new_data = zeros(FT, Nx, Ny, Nz)
        new_data .= (data[:, 1:end-1, :] .+ data[:, 2:end, :]) ./ 2
        new_data
    else
        Array{FT}(data)
    end

    FT_NaN = convert(FT, NaN)
    if metadata.name == :temperature && dataset_temperature_units(metadata) isa Kelvin
        data[data .!= FT_NaN] .-= FT(273.15) # convert to Celsius
    end

    set!(field, data)
end
=#