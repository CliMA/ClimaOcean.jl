#####
##### Zonal-mean sections and the equatorial undercurrent
#####

zonal_sections(output) = cached(output, :zonal_sections) do
    atlas = world_ocean_atlas(output)
    initial_buoyancy = record_array(timeseries(output, "fields", "bo"), 1)
    (; temperature = zonal_mean(output, time_mean(output, "fields", "to")),
       salinity = zonal_mean(output, time_mean(output, "fields", "so")),
       buoyancy = zonal_mean(output, time_mean(output, "fields", "bo")),
       atlas_temperature = zonal_mean(output, atlas.T),
       atlas_salinity = zonal_mean(output, atlas.S),
       initial_buoyancy = zonal_mean(output, initial_buoyancy))
end

finite_or_zero(x) = ifelse(isfinite(x), x, zero(x))

"""
The zonal root-mean-square of the time-mean temperature and salinity minus the World Ocean Atlas, and the part of it
that cancels in the zonal mean of the difference.
"""
zonal_rms_bias_sections(output) = cached(output, :zonal_rms_bias_sections) do
    atlas = world_ocean_atlas(output)
    sections = zonal_sections(output)
    rms(data, reference) = sqrt.(zonal_mean(output, finite_or_zero.(data .- reference) .^ 2))
    temperature = rms(time_mean(output, "fields", "to"), atlas.T)
    salinity = rms(time_mean(output, "fields", "so"), atlas.S)
    (; temperature, salinity,
       temperature_cancellation = temperature .- abs.(sections.temperature .- sections.atlas_temperature),
       salinity_cancellation = salinity .- abs.(sections.salinity .- sections.atlas_salinity))
end

zonal_mixed_layer_depth(output) = cached(output, :zonal_mixed_layer_depth) do
    monthly = abs.(cat(monthly_means(output, "surface", "mlotst")...; dims = 3))
    (; minimum = zonal_surface_mean(output, minimum(monthly; dims = 3)),
       maximum = zonal_surface_mean(output, maximum(monthly; dims = 3)))
end

function zonal_boyer_montegut_mixed_layer_depth()
    climatology = boyer_montegut_extrema()
    isnothing(climatology) && return nothing
    zonal(data) = [finite_mean(data[:, j]) for j in eachindex(map_latitudes)]
    return (; minimum = zonal(climatology.minimum), maximum = zonal(climatology.maximum))
end

#####
##### The equatorial undercurrent along 140ᵒW and along the equator
#####

longitude_distance(λ₁, λ₂) = abs(mod(λ₁ - λ₂ + 180, 360) - 180)

"""
The time-mean eastward velocity of `output` along 140ᵒW between 8ᵒS and 10ᵒN, and along the equator between 143ᵒE
and 95ᵒW, in the upper 400 m.
"""
equatorial_sections(output) = cached(output, :equatorial_sections) do
    u, v = centered_velocities(output, time_mean(output, "fields", "uo"), time_mean(output, "fields", "vo"))
    cosθ, sinθ = grid_rotation(output)
    uᵉ = u .* cosθ .- v .* sinθ
    λ, φ = center_longitudes(output), center_latitudes(output)
    Nx, Ny = size(λ)
    mask = ocean_mask(output)

    meridional = [(argmin(i -> longitude_distance(λ[i, j], 220), 1:Nx), j) for j in 1:Ny]
    meridional = filter(((i, j),) -> -8 ≤ φ[i, j] ≤ 10, meridional)

    pacific_column = argmin(i -> longitude_distance(λ[i, Ny ÷ 2], 180), 1:Nx)
    equator = argmin(j -> abs(φ[pacific_column, j]), 1:Ny)
    zonal = filter(i -> 143 ≤ λ[i, equator] ≤ 265, 1:Nx)

    section(cells) = permutedims(hcat((masked(uᵉ[i, j, :], mask[i, j, :]) for (i, j) in cells)...))

    (; latitudes = [φ[i, j] for (i, j) in meridional], along_140W = section(meridional),
       longitudes = λ[zonal, equator], along_equator = section([(i, equator) for i in zonal]),
       depths = depths(output))
end
