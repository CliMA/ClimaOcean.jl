#####
##### Sea ice concentration, extent, area and volume
#####

const extent_threshold = 0.15

sea_ice_concentration_maps(output) = cached(output, :sea_ice_concentration_maps) do
    monthly = monthly_means(output, "surface", "siconc")
    (; march = regrid_to_map(output, monthly[3]), september = regrid_to_map(output, monthly[9]))
end

function sea_ice_concentration_bias_maps(output)
    observed = hadisst_sea_ice_concentration()
    isnothing(observed) && return nothing
    model = sea_ice_concentration_maps(output)
    return (; march = model.march .- observed.march, september = model.september .- observed.september,
              observed_march = observed.march, observed_september = observed.september)
end

cell_areas(output) = cached(output, :cell_areas) do
    grid = model_grid(output)
    Nx, Ny, Nz = size(grid)
    [Azᶜᶜᶜ(i, j, Nz, grid) for i in 1:Nx, j in 1:Ny] .* surface_mask(output)
end

"""
The sea ice area, extent and volume of each hemisphere at every record of the surface output.
"""
sea_ice_integrals(output) = cached(output, :sea_ice_integrals) do
    ℵ = timeseries(output, "surface", "siconc")
    h = timeseries(output, "surface", "sithick")
    A = cell_areas(output)
    north = center_latitudes(output) .> 0
    south = .!north
    Nt = length(ℵ.times)

    integrals = (; arctic_area = zeros(Nt), antarctic_area = zeros(Nt),
                   arctic_extent = zeros(Nt), antarctic_extent = zeros(Nt),
                   arctic_volume = zeros(Nt), antarctic_volume = zeros(Nt))

    for n in 1:Nt
        ℵⁿ = dropdims(record_array(ℵ, n); dims = 3)
        hⁿ = dropdims(record_array(h, n); dims = 3)
        for (hemisphere, prefix) in ((north, :arctic), (south, :antarctic))
            getproperty(integrals, Symbol(prefix, :_area))[n] = sum(ℵⁿ .* A .* hemisphere)
            getproperty(integrals, Symbol(prefix, :_extent))[n] = sum((ℵⁿ .> extent_threshold) .* A .* hemisphere)
            getproperty(integrals, Symbol(prefix, :_volume))[n] = sum(hⁿ .* ℵⁿ .* A .* hemisphere)
        end
    end

    merge(integrals, (; years = in_years(ℵ.times), times = ℵ.times))
end

# The monthly climatology of a sea ice integral over the averaging window
function monthly_climatology(output, name)
    integrals = sea_ice_integrals(output)
    records = window_records(output, integrals.times)
    months = month.(record_dates(output, integrals.times[records]))
    values = getproperty(integrals, name)[records]
    return [mean(values[months .== m]) for m in 1:12]
end
