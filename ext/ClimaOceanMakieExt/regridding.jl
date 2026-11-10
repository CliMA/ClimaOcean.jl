#####
##### Conservative regridding onto the 1ᵒ grid every map and zonal mean is drawn on
#####

const map_grid = LatitudeLongitudeGrid(CPU(); size = (360, 180, 1), longitude = (0, 360), latitude = (-90, 90), z = (0, 1))
const map_longitudes = collect(λnodes(map_grid, Center()))
const map_latitudes = collect(φnodes(map_grid, Center()))

# Building a regridder takes minutes at high resolution, so runs on the same grid share one
const regridders = Dict{UInt, Any}()

regridder(output) = cached(output, :regridder) do
    key = hash((size(model_grid(output)), surface_mask(output)))
    get!(regridders, key) do
        @info "Building the regridder from the grid of $(output.label) to the 1ᵒ map grid..."
        ConservativeRegridding.Regridder(Field{Center, Center, Nothing}(map_grid),
                                         Field{Center, Center, Nothing}(model_grid(output)))
    end
end

# The mask-weighted regridding of one level, and the ocean area covering each map cell
function regrid_level(output, data, mask)
    weighted = zeros(length(map_longitudes) * length(map_latitudes))
    covered = zeros(length(weighted))
    valid = mask .& isfinite.(data)
    ConservativeRegridding.regrid!(weighted, regridder(output), vec(ifelse.(valid, data, 0.0)))
    ConservativeRegridding.regrid!(covered, regridder(output), vec(Float64.(valid)))
    shape = (length(map_longitudes), length(map_latitudes))
    return reshape(weighted, shape), reshape(covered, shape)
end

"""
The surface `data` of `output`, regridded to the 1ᵒ map grid; `NaN` where no ocean covers a map cell.
"""
function regrid_to_map(output, data, mask = surface_mask(output))
    weighted, covered = regrid_level(output, reshape(data, size(mask)), mask)
    return ifelse.(covered .> 0, weighted ./ covered, NaN)
end

"""
The zonal mean of the three-dimensional `data` of `output` on the latitudes of the map grid, as `(latitude, depth)`.
"""
function zonal_mean(output, data, mask = ocean_mask(output))
    Nz = size(data, 3)
    zonal = fill(NaN, length(map_latitudes), Nz)
    areas = reshape(regridder(output).dst_areas, length(map_longitudes), length(map_latitudes))

    for k in 1:Nz
        weighted, covered = regrid_level(output, view(data, :, :, k), view(mask, :, :, k))
        ocean_area = vec(sum(covered .* areas; dims = 1))
        content = vec(sum(weighted .* areas; dims = 1))
        zonal[:, k] .= ifelse.(ocean_area .> 0, content ./ ocean_area, NaN)
    end

    return zonal
end

zonal_surface_mean(output, data) =
    vec(zonal_mean(output, reshape(data, size(surface_mask(output))..., 1), ocean_mask(output)[:, :, end:end]))
