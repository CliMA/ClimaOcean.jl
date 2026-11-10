#####
##### Global integrals and horizontal means over the whole simulation
#####

const sea_ice_density = 900
const snow_density = 330
const sea_ice_salinity = 4

scalar_series(output, name) = cached(output, Symbol(:series, '/', name)) do
    fts = timeseries(output, "averages", name)
    [only(record_array(fts, n)) for n in eachindex(fts.times)]
end

series_years(output) = cached(() -> in_years(timeseries(output, "averages", "voco").times), output, :series_years)

ocean_volume(output) = scalar_series(output, "voco")
global_mean_temperature(output) = scalar_series(output, "hoco") ./ ocean_volume(output)
global_mean_salinity(output) = scalar_series(output, "soco") ./ ocean_volume(output)
global_mean_sea_level(output) = scalar_series(output, "zosga")

ocean_salt_content(output) = reference_density / 1000 .* scalar_series(output, "soco")
sea_ice_salt_content(output) = sea_ice_density * sea_ice_salinity / 1000 .* scalar_series(output, "sivol")
ocean_heat_content(output) = reference_density * heat_capacity .* scalar_series(output, "hoco")

total_water(output) = ocean_volume(output) .+
                      sea_ice_density / reference_density .* scalar_series(output, "sivol") .+
                      snow_density / reference_density .* scalar_series(output, "snvol")

profile_series(output, name) = cached(output, Symbol(:profiles, '/', name)) do
    fts = timeseries(output, "averages", name)
    permutedims(hcat((vec(record_array(fts, n)) for n in eachindex(fts.times))...))
end

function mean_profile(output, name)
    profiles = profile_series(output, name)
    records = window_records(output, timeseries(output, "averages", name).times)
    return vec(mean(profiles[records, :]; dims = 1))
end

profile_drift(output, name) = profile_series(output, name) .- profile_series(output, name)[1:1, :]

cell_volumes(output) = cached(output, :cell_volumes) do
    grid = model_grid(output)
    Nx, Ny, Nz = size(grid)
    [Vᶜᶜᶜ(i, j, k, grid) for i in 1:Nx, j in 1:Ny, k in 1:Nz] .* ocean_mask(output)
end

"""
The global-mean kinetic energy per unit mass of every record of the three-dimensional output.
"""
kinetic_energy_series(output) = cached(output, :kinetic_energy_series) do
    u² = timeseries(output, "fields", "uosq")
    v² = timeseries(output, "fields", "vosq")
    V = cell_volumes(output)
    energy = [sum((record_array(u², n) .+ record_array(v², n)) .* V) / 2sum(V) for n in eachindex(u².times)]
    (; years = in_years(u².times), energy)
end

"""
The horizontal root-mean-square difference between the temperature and salinity of every record of the
three-dimensional output and the World Ocean Atlas, at every level.
"""
rms_bias_series(output) = cached(output, :rms_bias_series) do
    T = timeseries(output, "fields", "to")
    S = timeseries(output, "fields", "so")
    atlas = world_ocean_atlas(output)
    V = cell_volumes(output) .* isfinite.(atlas.T) .* isfinite.(atlas.S)
    level_volume = vec(sum(V; dims = (1, 2)))
    finite(x) = ifelse(isfinite(x), x, 0.0)

    rms(data, reference) = sqrt.(vec(sum(finite.(data .- reference) .^ 2 .* V; dims = (1, 2))) ./ level_volume)

    Nt = length(T.times)
    temperature = permutedims(hcat((rms(record_array(T, n), atlas.T) for n in 1:Nt)...))
    salinity = permutedims(hcat((rms(record_array(S, n), atlas.S) for n in 1:Nt)...))

    (; years = in_years(T.times), temperature, salinity)
end
