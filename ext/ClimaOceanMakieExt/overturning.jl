#####
##### Overturning and barotropic streamfunctions
#####

const rapid_latitude = 26.5
const rapid_index_minimum_depth = 500

atlantic_mask(output) = cached(output, :atlantic_mask) do
    dropdims(Array(interior(atlantic_ocean_basin(model_grid(output)).mask)); dims = 3) .> 0
end

# The latitude of every row of meridional velocity faces, averaged along the row
overturning_latitudes(output, Ny) = cached(output, :overturning_latitudes) do
    grid = model_grid(output)
    Nx = size(grid, 1)
    [mean(φnode(i, j, 1, grid, Center(), Face(), Center()) for i in 1:Nx) for j in 1:Ny]
end

"""
The Atlantic overturning streamfunction [Sv] of the meridional volume transport `transport`, integrated from the
bottom, as `(row, depth)`.
"""
function atlantic_streamfunction(output, transport)
    mask = atlantic_mask(output)
    Nx, Ny = size(mask)
    zonal_transport = dropdims(sum(transport[:, 1:Ny, :] .* mask; dims = 1); dims = 1)
    return - cumsum(zonal_transport; dims = 2) ./ 1e6
end

has_bolus_transport(output) = has_timeseries(output, "fields", "vvolgm")

atlantic_overturning_streamfunction(output) = cached(output, :atlantic_overturning) do
    ψ = atlantic_streamfunction(output, time_mean(output, "fields", "vvol"))
    has_bolus_transport(output) && (ψ .+= atlantic_streamfunction(output, time_mean(output, "fields", "vvolgm")))
    (; ψ, latitudes = overturning_latitudes(output, size(ψ, 1)), depths = depths(output))
end

rapid_row(output) = argmin(abs.(atlantic_overturning_streamfunction(output).latitudes .- rapid_latitude))

"""
The maximum of the Atlantic overturning streamfunction at 26.5ᵒN below 500 m at every record of the
three-dimensional output.
"""
rapid_index_series(output) = cached(output, :rapid_index_series) do
    v = timeseries(output, "fields", "vvol")
    v★ = has_bolus_transport(output) ? timeseries(output, "fields", "vvolgm") : nothing
    j = rapid_row(output)
    deep = depths(output) .< -rapid_index_minimum_depth

    ψ = map(eachindex(v.times)) do n
        transport = record_array(v, n)
        isnothing(v★) || (transport .+= record_array(v★, n))
        maximum(atlantic_streamfunction(output, transport)[j, deep])
    end

    (; years = in_years(v.times), ψ)
end

barotropic_streamfunction_map(output) = cached(output, :barotropic_streamfunction) do
    depth_integrated = dropdims(sum(time_mean(output, "fields", "uvol"); dims = 3); dims = 3)
    Ψ = - cumsum(depth_integrated; dims = 2) ./ 1e6
    λ, φ = center_longitudes(output), center_latitudes(output)
    land = findall(.!surface_mask(output))
    # Referenced to a land point of South America, so that the gyres circulate around zero
    reference = argmin(I -> longitude_distance(λ[I], 300)^2 + (φ[I] + 10)^2, land)
    regrid_to_map(output, Ψ .- Ψ[reference])
end
