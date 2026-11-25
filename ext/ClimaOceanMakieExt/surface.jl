#####
##### Surface climatologies and their biases, on the map grid
#####

const reference_density = 1026
const heat_capacity = 3991.86795711963

surface_level(data) = data[:, :, end:end]

sea_surface_temperature_bias_map(output) = cached(output, :sea_surface_temperature_bias) do
    regrid_to_map(output, time_mean(output, "surface", "tos") .- surface_level(world_ocean_atlas(output).T))
end

sea_surface_salinity_bias_map(output) = cached(output, :sea_surface_salinity_bias) do
    regrid_to_map(output, time_mean(output, "surface", "sos") .- surface_level(world_ocean_atlas(output).S))
end

sea_surface_height_map(output) = cached(() -> regrid_to_map(output, time_mean(output, "surface", "zos")), output, :sea_surface_height)

area_weighted_mean(data) = let weights = cosd.(map_latitudes') .* isfinite.(data)
    sum(ifelse.(isfinite.(data), data, 0) .* weights) / sum(weights)
end

demeaned(data) = data .- area_weighted_mean(data)

function sea_surface_height_difference_map(output)
    ecco = ecco_sea_surface_height()
    isnothing(ecco) && return nothing
    model = sea_surface_height_map(output)
    return demeaned(model) .- demeaned(ifelse.(isfinite.(model), ecco, NaN))
end

sea_surface_height_variability_map(output) = cached(output, :sea_surface_height_variability) do
    η = time_mean(output, "surface", "zos")
    η² = time_mean(output, "surface", "zossq")
    regrid_to_map(output, sqrt.(max.(η² .- η .^ 2, 0)))
end

# Summer minimum and winter maximum of the monthly mean mixed-layer depth
mixed_layer_depth_extrema_maps(output) = cached(output, :mixed_layer_depth_extrema) do
    monthly = abs.(cat(monthly_means(output, "surface", "mlotst")...; dims = 3))
    (; minimum = regrid_to_map(output, minimum(monthly; dims = 3)),
       maximum = regrid_to_map(output, maximum(monthly; dims = 3)))
end

boyer_montegut_extrema() = let monthly = boyer_montegut_mixed_layer_depth()
    isnothing(monthly) ? nothing : (; minimum = minimum(cat(monthly...; dims = 3); dims = 3)[:, :, 1],
                                       maximum = maximum(cat(monthly...; dims = 3); dims = 3)[:, :, 1])
end

# The fluxes are written positive out of the ocean
net_heat_flux_map(output) = cached(output, :net_heat_flux) do
    regrid_to_map(output, - reference_density * heat_capacity .* time_mean(output, "surface", "hfds"))
end

net_freshwater_flux_map(output) = cached(output, :net_freshwater_flux) do
    regrid_to_map(output, reference_density .* time_mean(output, "surface", "wfo") ./ time_mean(output, "surface", "sos"))
end

wind_stress_maps(output) = cached(output, :wind_stress) do
    τx, τy = centered_velocities(output, - reference_density .* time_mean(output, "surface", "tauuo"),
                                         - reference_density .* time_mean(output, "surface", "tauvo"))
    τᵉ, τⁿ = eastward_northward(output, τx, τy)
    (; eastward = regrid_to_map(output, τᵉ), northward = regrid_to_map(output, τⁿ))
end

near_surface_current_maps(output) = cached(output, :near_surface_currents) do
    u, v = centered_velocities(output, time_mean(output, "surface", "uos"), time_mean(output, "surface", "vos"))
    uᵉ, vⁿ = eastward_northward(output, u, v)
    eastward = regrid_to_map(output, uᵉ)
    northward = regrid_to_map(output, vⁿ)
    (; eastward, northward, speed = hypot.(eastward, northward))
end
