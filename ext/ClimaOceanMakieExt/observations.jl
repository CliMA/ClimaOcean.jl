#####
##### Observations, on the map grid or as tables
#####

const observations = Dict{Any, Any}()

# Observations are a reference: a figure is drawn without them when they cannot be obtained
function observation(f, key)
    return get!(observations, key) do
        try
            f()
        catch error
            @warn "$(first(key)) is not available: the figures are drawn without it" exception = (error, catch_backtrace())
            nothing
        end
    end
end

observations_directory() = mkpath(Observations.observations_cache)

# Reductions of long records are computed once and kept next to the downloads
function saved(f, name)
    path = joinpath(observations_directory(), name * ".jld2")
    isfile(path) && return jldopen(file -> file["value"], path)
    value = f()
    jldopen(file -> (file["value"] = value), path, "w")
    return value
end

function on_map_grid(metadatum)
    field = Field{Center, Center, Nothing}(map_grid)
    set!(field, metadatum; inpainting = nothing)
    return dropdims(Array(interior(field)); dims = 3)
end

monthly_maps(name, dataset) = [on_map_grid(Metadatum(name; dataset, date)) for date in all_dates(dataset, name)]

#####
##### World Ocean Atlas on the grid of a run, in Conservative Temperature and Absolute Salinity
#####

world_ocean_atlas(output) = cached(output, :world_ocean_atlas) do
    grid = model_grid(output)
    T = CenterField(grid)
    S = CenterField(grid)
    set!(T, Metadatum(:temperature; dataset = WOAAnnual()))
    set!(S, Metadatum(:salinity; dataset = WOAAnnual()))
    convert_to_teos10!(T, S)
    mask = ocean_mask(output)
    (; T = masked(Array(interior(T)), mask), S = masked(Array(interior(S)), mask))
end

#####
##### Gridded climatologies
#####

ecco_sea_surface_height() = observation(("ECCO sea surface height",)) do
    saved("ecco_sea_surface_height_1992_2012") do
        dates = DateTime(1992, 1, 1):Dates.Month(1):DateTime(2012, 12, 1)
        sum(on_map_grid(Metadatum(:free_surface; dataset = ECCO4Monthly(), date)) for date in dates) ./ length(dates)
    end
end

boyer_montegut_mixed_layer_depth() = observation(("de Boyer Montégut mixed-layer depth",)) do
    monthly_maps(:mixed_layer_depth, BoyerMontegutMonthly())
end

ncep_wind_stress() = observation(("NCEP wind stress",)) do
    (sum(monthly_maps(:zonal_wind_stress, NCEPMonthly())) ./ 12,
     sum(monthly_maps(:meridional_wind_stress, NCEPMonthly())) ./ 12)
end

const hadisst_years = (1979, 2007)

# March and September of the HadISST sea ice concentration averaged over `hadisst_years`
hadisst_sea_ice_concentration() = observation(("HadISST sea ice concentration",)) do
    saved("hadisst_sea_ice_concentration_$(hadisst_years[1])_$(hadisst_years[2])") do
        climatology(m) = sum(on_map_grid(Metadatum(:sea_ice_concentration; dataset = HadISSTMonthly(), date = DateTime(y, m, 1)))
                             for y in hadisst_years[1]:hadisst_years[2]) ./ (hadisst_years[2] - hadisst_years[1] + 1)
        (; march = climatology(3), september = climatology(9))
    end
end

#####
##### Sea ice extent and area of NSIDC, and Arctic sea ice volume of PIOMAS
#####

function downloaded(url)
    path = joinpath(observations_directory(), basename(url))
    isfile(path) || download_with_retries(url, path)
    return path
end

function monthly_mean_and_spread(years, values, window)
    keep = window[1] .≤ years .≤ window[2]
    any(keep) || error("no observation between $(window[1]) and $(window[2])")
    return mean(values[keep]), count(keep) > 1 ? std(values[keep]) : 0.0, extrema(years[keep])
end

nsidc(hemisphere, window) = observation(("NSIDC sea ice index", hemisphere, window)) do
    letter = hemisphere == :north ? "N" : "S"
    records = map(1:12) do m
        url = "https://noaadata.apps.nsidc.org/NOAA/G02135/$hemisphere/monthly/data/$(letter)_$(lpad(m, 2, '0'))_extent_v4.0.csv"
        rows = [split(line, ',') for line in readlines(downloaded(url))]
        rows = [row for row in rows if length(row) ≥ 6 && !isnothing(tryparse(Int, strip(row[1])))]
        years = [parse(Int, strip(row[1])) for row in rows]
        extent = [parse(Float64, strip(row[5])) for row in rows]
        area = [parse(Float64, strip(row[6])) for row in rows]
        valid = (extent .> 0) .& (area .> 0)
        (years[valid], extent[valid], area[valid])
    end

    extent = [monthly_mean_and_spread(years, extent, window) for (years, extent, _) in records]
    area = [monthly_mean_and_spread(years, area, window) for (years, _, area) in records]

    (; extent = first.(extent), extent_spread = getindex.(extent, 2),
       area = first.(area), area_spread = getindex.(area, 2), years = last(first(extent)))
end

const piomas_url = "https://psc.apl.uw.edu/wordpress/wp-content/uploads/schweiger/ice_volume/PIOMAS.monthly.Current.v2.1.csv"

piomas(window) = observation(("PIOMAS Arctic sea ice volume", window)) do
    table = readdlm(downloaded(piomas_url), ','; skipstart = 1)
    years = Int.(table[:, 1])
    volume = Float64.(table[:, 2:13])
    statistics = [monthly_mean_and_spread(years[volume[:, m] .> 0], volume[volume[:, m] .> 0, m], window) for m in 1:12]
    (; volume = first.(statistics), volume_spread = getindex.(statistics, 2), years = last(first(statistics)))
end

observation_years_label(years) = string(years[1], "–", years[2])

#####
##### The Atlantic overturning at 26.5ᵒN measured by RAPID
#####

const rapid_profile_url = "https://rapid.ac.uk/sites/default/files/rapid_data/moc_vertical.nc"
const rapid_index_url = "https://rapid.ac.uk/sites/default/files/rapid_data/moc_transports.nc"

finite_mean(values) = let finite = filter(isfinite, values)
    isempty(finite) ? NaN : mean(finite)
end

finite_spread(values) = let finite = filter(isfinite, values)
    length(finite) < 2 ? NaN : std(finite)
end

rapid_profile() = observation(("RAPID overturning profile",)) do
    NCDataset(downloaded(rapid_profile_url)) do dataset
        depth = Float64.(dataset["depth"][:])
        ψ = Float64.(coalesce.(dataset["stream_function_mar"][:, :], NaN))
        ψ = size(ψ, 1) == length(depth) ? ψ : permutedims(ψ)
        (; z = -depth, ψ = [finite_mean(ψ[k, :]) for k in eachindex(depth)],
           spread = [finite_spread(ψ[k, :]) for k in eachindex(depth)])
    end
end

rapid_index() = observation(("RAPID overturning index",)) do
    NCDataset(downloaded(rapid_index_url)) do dataset
        dates = dataset["time"][:]
        ψ = Float64.(coalesce.(dataset["moc_mar_hc10"][:], NaN))
        (; year = year.(dates) .+ (Dates.dayofyear.(dates) .- 0.5) ./ 365.25, ψ)
    end
end
