#####
##### The OMIP figures
#####

"""
Rows of maps, one column per run: each row draws `data(output)` with `draw`, or leaves the panel empty where `data`
returns `nothing`. An `observation` row is drawn once, in the first column.
"""
function maps_figure(outputs, rows)
    columns = length(outputs)
    return omip_figure(columns, length(rows)) do figure
        for (r, row) in enumerate(rows)
            plot = nothing
            for (c, output) in enumerate(row.observation ? outputs[1:1] : outputs)
                data = row.data(output)
                isnothing(data) && continue
                title = row.observation ? row.title : string(output.label, ": ", row.title)
                plot = row.draw(panel(figure, r, c), data; title, colormap = row.colormap, colorrange = row.colorrange)
            end
            isnothing(plot) || row_colorbar!(figure, r, columns, plot, row.label)
        end
    end
end

map_row(title, data, colormap, colorrange, label; draw = global_map!, observation = false) =
    (; title, data, colormap, colorrange, label, draw, observation)

ClimaOcean.sea_surface_temperature_bias(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("SST − WOA", sea_surface_temperature_bias_map, :balance, (-3, 3), "ᵒC")])

ClimaOcean.sea_surface_salinity_bias(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("SSS − WOA", sea_surface_salinity_bias_map, :balance, (-2, 2), "g kg⁻¹")])

ClimaOcean.sea_surface_height(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("time-mean SSH", sea_surface_height_map, :balance, (-1.5, 1.5), "m"),
                          map_row("SSH − ECCO 1992–2012, both demeaned", sea_surface_height_difference_map,
                                  :balance, (-0.5, 0.5), "m")])

ClimaOcean.sea_surface_height_variability(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("SSH variability", sea_surface_height_variability_map, :magma, (0, 0.2), "m")])

ClimaOcean.surface_fluxes(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("net heat flux into the ocean", net_heat_flux_map, :balance, (-200, 200), "W m⁻²"),
                          map_row("net freshwater flux into the ocean", net_freshwater_flux_map, :balance,
                                  (-8e-5, 8e-5), "kg m⁻² s⁻¹")])

function ncep_difference(output, component)
    observed = ncep_wind_stress()
    isnothing(observed) && return nothing
    return getproperty(wind_stress_maps(output), component) .- observed[component == :eastward ? 1 : 2]
end

ClimaOcean.wind_stress(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("eastward wind stress", o -> wind_stress_maps(o).eastward, :balance, (-0.2, 0.2), "N m⁻²"),
                          map_row("northward wind stress", o -> wind_stress_maps(o).northward, :balance, (-0.2, 0.2), "N m⁻²"),
                          map_row("eastward stress − NCEP", o -> ncep_difference(o, :eastward), :balance, (-0.1, 0.1), "N m⁻²"),
                          map_row("northward stress − NCEP", o -> ncep_difference(o, :northward), :balance, (-0.1, 0.1), "N m⁻²")])

ClimaOcean.near_surface_currents(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("near-surface speed", o -> near_surface_current_maps(o).speed, :speed, (0, 0.5), "m s⁻¹"),
                          map_row("eastward velocity", o -> near_surface_current_maps(o).eastward, :balance, (-0.5, 0.5), "m s⁻¹"),
                          map_row("northward velocity", o -> near_surface_current_maps(o).northward, :balance, (-0.3, 0.3), "m s⁻¹")])

const mixed_layer_colormap = Reverse(:deep)

function ClimaOcean.mixed_layer_depth(outputs::OMIPOutput...)
    rows = [map_row("summer minimum MLD", o -> mixed_layer_depth_extrema_maps(o).minimum, mixed_layer_colormap, (0, 80), "m"),
            map_row("winter maximum MLD", o -> mixed_layer_depth_extrema_maps(o).maximum, mixed_layer_colormap, (0, 600), "m")]

    observed = boyer_montegut_extrema()
    isnothing(observed) ||
        append!(rows, [map_row("de Boyer Montégut et al. (2004): summer minimum", _ -> observed.minimum,
                               mixed_layer_colormap, (0, 80), "m"; observation = true),
                       map_row("de Boyer Montégut et al. (2004): winter maximum", _ -> observed.maximum,
                               mixed_layer_colormap, (0, 600), "m"; observation = true)])

    return maps_figure(outputs, rows)
end

const polar_panels = (("March, Arctic", :march, :north), ("March, Antarctic", :march, :south),
                      ("September, Arctic", :september, :north), ("September, Antarctic", :september, :south))

polar_map_row(title, data, hemisphere, colormap, colorrange, label; edge = nothing) =
    (; title, data, colormap, colorrange, label, hemisphere, edge)

function polar_figure(outputs, rows)
    columns = length(outputs)
    return omip_figure(columns, length(rows); column_width = 420, row_height = 420) do figure
        for (r, row) in enumerate(rows)
            plot = nothing
            for (c, output) in enumerate(outputs)
                data = row.data(output)
                isnothing(data) && continue
                plot = polar_map!(panel(figure, r, c), data; hemisphere = row.hemisphere, edge = row.edge(output),
                                  title = string(output.label, ": ", row.title),
                                  colormap = row.colormap, colorrange = row.colorrange)
            end
            isnothing(plot) || row_colorbar!(figure, r, columns, plot, row.label)
        end
    end
end

ClimaOcean.sea_ice_concentration(outputs::OMIPOutput...) =
    polar_figure(outputs, [(; title, hemisphere, colormap = :ice, colorrange = (0, 1), label = "concentration",
                              data = o -> getproperty(sea_ice_concentration_maps(o), season), edge = o -> nothing)
                           for (title, season, hemisphere) in polar_panels])

function ClimaOcean.sea_ice_concentration_bias(outputs::OMIPOutput...)
    observed_edge(season) = output -> let bias = sea_ice_concentration_bias_maps(output)
        isnothing(bias) ? nothing : getproperty(bias, Symbol(:observed_, season))
    end
    bias(season) = output -> let bias = sea_ice_concentration_bias_maps(output)
        isnothing(bias) ? nothing : getproperty(bias, season)
    end
    return polar_figure(outputs, [(; title = "$title − HadISST", hemisphere, colormap = Reverse(:RdBu),
                                     colorrange = (-0.6, 0.6), label = "concentration difference, observed 15% edge in black",
                                     data = bias(season), edge = observed_edge(season))
                                  for (title, season, hemisphere) in polar_panels])
end

#####
##### Climatologies and time series
#####

function lines_figure(f, panels; legend = true)
    return omip_figure(length(panels) + legend, 1; column_width = 520, row_height = 420) do figure
        axes = [f_axis(figure[1, p]) for (p, f_axis) in enumerate(panels)]
        f(axes)
        legend && legend!(figure[1, length(panels) + 1], first(axes))
    end
end

const square_meters_to_million_square_kilometers = 1e-12
const cubic_meters_to_thousand_cubic_kilometers = 1e-12

function sea_ice_climatology_figure(outputs, quantity, ylabel)
    first_output = first(outputs)
    window = averaged_forcing_years(first_output)
    panels = [position -> climatology_axis(position; title = "Arctic", ylabel),
              position -> climatology_axis(position; title = "Antarctic", ylabel)]
    return lines_figure(panels) do axes
        for (axis, hemisphere, prefix) in zip(axes, (:north, :south), (:arctic, :antarctic))
            observed = nsidc(hemisphere, window)
            isnothing(observed) ||
                observation_lines!(axis, 1:12, getproperty(observed, quantity), getproperty(observed, Symbol(quantity, :_spread)),
                                   "NSIDC " * observation_years_label(observed.years))
            for (n, output) in enumerate(outputs)
                run_lines!(axis, 1:12, monthly_climatology(output, Symbol(prefix, '_', quantity)) .*
                                        square_meters_to_million_square_kilometers, n, output.label)
            end
        end
    end
end

ClimaOcean.sea_ice_extent(outputs::OMIPOutput...) = sea_ice_climatology_figure(outputs, :extent, "sea ice extent (10⁶ km²)")
ClimaOcean.sea_ice_area(outputs::OMIPOutput...) = sea_ice_climatology_figure(outputs, :area, "sea ice area (10⁶ km²)")

function ClimaOcean.arctic_sea_ice_volume(outputs::OMIPOutput...)
    observed = piomas(averaged_forcing_years(first(outputs)))
    return lines_figure([position -> climatology_axis(position; title = "Arctic sea ice volume", ylabel = "volume (10³ km³)")]) do (axis,)
        isnothing(observed) || observation_lines!(axis, 1:12, observed.volume, observed.volume_spread,
                                                  "PIOMAS " * observation_years_label(observed.years))
        for (n, output) in enumerate(outputs)
            run_lines!(axis, 1:12, monthly_climatology(output, :arctic_volume) .* cubic_meters_to_thousand_cubic_kilometers,
                       n, output.label)
        end
    end
end

function series_figure(outputs, series; titles, ylabels)
    panels = [position -> time_axis(position; title, ylabel) for (title, ylabel) in zip(titles, ylabels)]
    return lines_figure(panels) do axes
        for (axis, s) in zip(axes, series), (n, output) in enumerate(outputs)
            years, values = s(output)
            run_lines!(axis, years, values, n, output.label)
        end
    end
end

ClimaOcean.sea_ice_area_time_series(outputs::OMIPOutput...) =
    series_figure(outputs, [o -> (sea_ice_integrals(o).years, sea_ice_integrals(o).arctic_area .* square_meters_to_million_square_kilometers),
                            o -> (sea_ice_integrals(o).years, sea_ice_integrals(o).antarctic_area .* square_meters_to_million_square_kilometers)];
                  titles = ("Arctic sea ice area", "Antarctic sea ice area"), ylabels = ("area (10⁶ km²)", "area (10⁶ km²)"))

ClimaOcean.arctic_sea_ice_volume_time_series(outputs::OMIPOutput...) =
    series_figure(outputs, [o -> (sea_ice_integrals(o).years, sea_ice_integrals(o).arctic_volume .* cubic_meters_to_thousand_cubic_kilometers)];
                  titles = ("Arctic sea ice volume",), ylabels = ("volume (10³ km³)",))

ClimaOcean.kinetic_energy(outputs::OMIPOutput...) =
    series_figure(outputs, [o -> (kinetic_energy_series(o).years, kinetic_energy_series(o).energy)];
                  titles = ("global-mean kinetic energy",), ylabels = ("kinetic energy (m² s⁻²)",))

drift(values) = values .- first(values)

ClimaOcean.global_drift(outputs::OMIPOutput...) =
    series_figure(outputs, [o -> (series_years(o), drift(global_mean_temperature(o))),
                            o -> (series_years(o), drift(global_mean_salinity(o)))];
                  titles = ("global-mean temperature drift", "global-mean salinity drift"),
                  ylabels = ("ΔT (ᵒC)", "ΔS (g kg⁻¹)"))

ClimaOcean.sea_level_drift(outputs::OMIPOutput...) =
    series_figure(outputs, [o -> (series_years(o), global_mean_sea_level(o))];
                  titles = ("global-mean sea surface height",), ylabels = ("⟨η⟩ (m)",))

function ClimaOcean.content_conservation(outputs::OMIPOutput...)
    series = [o -> (series_years(o), drift(ocean_salt_content(o)) ./ 1e12),
              o -> (series_years(o), drift(sea_ice_salt_content(o)) ./ 1e12),
              o -> (series_years(o), drift(ocean_salt_content(o) .+ sea_ice_salt_content(o)) ./ 1e12),
              o -> (series_years(o), drift(ocean_heat_content(o)) ./ 1e21),
              o -> (series_years(o), drift(ocean_volume(o)) ./ 1e9),
              o -> (series_years(o), drift(total_water(o)) ./ 1e9)]
    titles = ("ocean salt", "sea ice salt", "ocean and sea ice salt", "ocean heat", "ocean volume",
              "ocean, sea ice and snow water")
    ylabels = ("Δ salt (Gt)", "Δ salt (Gt)", "Δ salt (Gt)", "Δ heat (ZJ)", "Δ volume (km³)", "Δ volume (km³)")

    return omip_figure(4, 2; column_width = 480, row_height = 400) do figure
        axes = [time_axis(figure[fldmod1(p, 3)...]; title, ylabel) for (p, (title, ylabel)) in enumerate(zip(titles, ylabels))]
        for (axis, s) in zip(axes, series)
            hlines!(axis, 0; color = (:black, 0.3), linewidth = 1)
            for (n, output) in enumerate(outputs)
                years, values = s(output)
                run_lines!(axis, years, values, n, output.label)
            end
        end
        legend!(figure[1:2, 4], first(axes))
    end
end

function ClimaOcean.mean_profiles(outputs::OMIPOutput...)
    panels = [position -> Axis(position; title = "horizontal-mean temperature", xlabel = "T (ᵒC)", ylabel = "Depth (m)"),
              position -> Axis(position; title = "horizontal-mean salinity", xlabel = "S (g kg⁻¹)", ylabel = "Depth (m)")]
    return lines_figure(panels) do axes
        for (axis, name) in zip(axes, ("to_h", "so_h")), (n, output) in enumerate(outputs)
            run_lines!(axis, mean_profile(output, name), depths(output), n, output.label)
            ylims!(axis, -maximum_depth, 0)
        end
    end
end

function ClimaOcean.zonal_mean_mixed_layer_depth(outputs::OMIPOutput...)
    observed = zonal_boyer_montegut_mixed_layer_depth()
    panels = [position -> Axis(position; title = "zonal-mean summer minimum", xlabel = "Latitude", ylabel = "MLD (m)"),
              position -> Axis(position; title = "zonal-mean winter maximum", xlabel = "Latitude", ylabel = "MLD (m)")]
    return lines_figure(panels) do axes
        for (axis, extremum) in zip(axes, (:minimum, :maximum))
            isnothing(observed) || lines!(axis, map_latitudes, getproperty(observed, extremum);
                                          observation_line..., label = "de Boyer Montégut et al. (2004)")
            for (n, output) in enumerate(outputs)
                run_lines!(axis, map_latitudes, getproperty(zonal_mixed_layer_depth(output), extremum), n, output.label)
            end
        end
    end
end

#####
##### Sections
#####

function sections_figure(outputs, rows)
    columns = length(outputs)
    return omip_figure(columns, length(rows); column_width = 560, row_height = 330) do figure
        for (r, row) in enumerate(rows)
            plot = nothing
            for (c, output) in enumerate(outputs)
                x, z, data, overlay = row.data(output)
                plot = section_panel!(panel(figure, r, c), x, z, data; title = string(output.label, ": ", row.title),
                                      levels = row.levels, colormap = row.colormap, overlay, xlabel = row.xlabel,
                                      depth_limit = row.depth_limit)
            end
            row_colorbar!(figure, r, columns, plot, row.label)
        end
    end
end

section_row(title, data, levels, colormap, label; xlabel = "Latitude", depth_limit = maximum_depth) =
    (; title, data, levels, colormap, label, xlabel, depth_limit)

zonal_section(output, name, overlay_name = nothing) =
    let sections = zonal_sections(output)
        (map_latitudes, depths(output), getproperty(sections, name),
         isnothing(overlay_name) ? nothing : getproperty(sections, overlay_name))
    end

ClimaOcean.zonal_mean_sections(outputs::OMIPOutput...) =
    sections_figure(outputs, [section_row("zonal-mean T, WOA dashed", o -> zonal_section(o, :temperature, :atlas_temperature),
                                          -2:2:30, :thermal, "ᵒC"),
                              section_row("zonal-mean S, WOA dashed", o -> zonal_section(o, :salinity, :atlas_salinity),
                                          33:0.25:37, :haline, "g kg⁻¹"),
                              section_row("zonal-mean b, initial dashed", o -> zonal_section(o, :buoyancy, :initial_buoyancy),
                                          range(-0.04, 0.02, length = 13), :balance, "m s⁻²")])

function zonal_difference(output, name, reference)
    sections = zonal_sections(output)
    return (map_latitudes, depths(output), getproperty(sections, name) .- getproperty(sections, reference), nothing)
end

function ClimaOcean.zonal_mean_drift(outputs::OMIPOutput...)
    largest_buoyancy_drift = maximum(o -> maximum(abs, filter(isfinite, zonal_difference(o, :buoyancy, :initial_buoyancy)[3])), outputs)
    return sections_figure(outputs, [section_row("zonal-mean T − WOA", o -> zonal_difference(o, :temperature, :atlas_temperature),
                                                 range(-3, 3, length = 13), :balance, "ᵒC"),
                                     section_row("zonal-mean S − WOA", o -> zonal_difference(o, :salinity, :atlas_salinity),
                                                 range(-0.75, 0.75, length = 13), :balance, "g kg⁻¹"),
                                     section_row("zonal-mean b − initial", o -> zonal_difference(o, :buoyancy, :initial_buoyancy),
                                                 range(-largest_buoyancy_drift, largest_buoyancy_drift, length = 13), :balance, "m s⁻²")])
end

rms_section(output, name) = (map_latitudes, depths(output), getproperty(zonal_rms_bias_sections(output), name), nothing)

ClimaOcean.zonal_rms_bias(outputs::OMIPOutput...) =
    sections_figure(outputs, [section_row("zonal RMS of T − WOA", o -> rms_section(o, :temperature), range(0, 3, length = 13), :thermal, "ᵒC"),
                              section_row("zonal RMS of S − WOA", o -> rms_section(o, :salinity), range(0, 0.75, length = 13), :haline, "g kg⁻¹"),
                              section_row("T bias cancelled in the zonal mean", o -> rms_section(o, :temperature_cancellation),
                                          range(0, 3, length = 13), :thermal, "ᵒC"),
                              section_row("S bias cancelled in the zonal mean", o -> rms_section(o, :salinity_cancellation),
                                          range(0, 0.75, length = 13), :haline, "g kg⁻¹")])

function ClimaOcean.atlantic_overturning(outputs::OMIPOutput...)
    overturning(output) = let ψ = atlantic_overturning_streamfunction(output)
        (ψ.latitudes, ψ.depths, ψ.ψ, nothing)
    end
    figure = sections_figure(outputs, [section_row("Atlantic overturning", overturning, -6:2:20, :balance, "Sv")])
    for content in contents(figure.layout)
        content isa Axis && xlims!(content, -35, 70)
    end
    return figure
end

function ClimaOcean.equatorial_undercurrent(outputs::OMIPOutput...)
    along_140W(output) = let s = equatorial_sections(output)
        (s.latitudes, s.depths, 100 .* s.along_140W, nothing)
    end
    along_equator(output) = let s = equatorial_sections(output)
        (s.longitudes, s.depths, 100 .* s.along_equator, nothing)
    end
    return sections_figure(outputs, [section_row("eastward velocity along 140ᵒW", along_140W, -100:10:100, :balance, "cm s⁻¹";
                                                 depth_limit = 400),
                                     section_row("eastward velocity along the equator", along_equator, -100:10:100, :balance,
                                                 "cm s⁻¹"; xlabel = "Longitude (ᵒE)", depth_limit = 400)])
end

function depth_time_figure(outputs, rows)
    columns = length(outputs)
    return omip_figure(columns, length(rows); column_width = 720, row_height = 520) do figure
        for (r, row) in enumerate(rows)
            plot = nothing
            for (c, output) in enumerate(outputs)
                years, data = row.data(output)
                plot = depth_time_panel!(panel(figure, r, c), years, depths(output), data;
                                         title = string(output.label, ": ", row.title), levels = row.levels,
                                         colormap = row.colormap)
            end
            row_colorbar!(figure, r, columns, plot, row.label)
        end
    end
end

ClimaOcean.drift_profiles(outputs::OMIPOutput...) =
    depth_time_figure(outputs, [(; title = "horizontal-mean ΔT", data = o -> (series_years(o), profile_drift(o, "to_h")),
                                   levels = range(-1.6, 1.6, length = 17), colormap = :balance, label = "ᵒC"),
                                (; title = "horizontal-mean ΔS", data = o -> (series_years(o), profile_drift(o, "so_h")),
                                   levels = range(-0.1, 0.1, length = 21), colormap = :balance, label = "g kg⁻¹")])

ClimaOcean.rms_bias_profiles(outputs::OMIPOutput...) =
    depth_time_figure(outputs, [(; title = "horizontal RMS of T − WOA", data = o -> (rms_bias_series(o).years, rms_bias_series(o).temperature),
                                   levels = range(0, 3, length = 13), colormap = :thermal, label = "ᵒC"),
                                (; title = "horizontal RMS of S − WOA", data = o -> (rms_bias_series(o).years, rms_bias_series(o).salinity),
                                   levels = range(0, 0.6, length = 13), colormap = :haline, label = "g kg⁻¹")])

#####
##### Overturning at 26.5ᵒN, streamfunctions and transports
#####

function six_month_means(years, values)
    half_years = floor.(2 .* years) ./ 2
    centers = sort(unique(half_years))
    return centers .+ 0.25, [finite_mean(values[half_years .== h]) for h in centers]
end

function ClimaOcean.atlantic_overturning_at_rapid(outputs::OMIPOutput...)
    profile, index = rapid_profile(), rapid_index()
    panels = [position -> Axis(position; title = "Atlantic overturning at 26.5ᵒN", xlabel = "ψ (Sv)", ylabel = "Depth (m)"),
              position -> Axis(position; title = "maximum below 500 m, six-month means", xlabel = "Year", ylabel = "ψ (Sv)")]

    return lines_figure(panels) do (profile_axis, index_axis)
        if !isnothing(profile)
            band!(profile_axis, Point2f.(profile.ψ .- profile.spread, profile.z), Point2f.(profile.ψ .+ profile.spread, profile.z);
                  observation_band...)
            lines!(profile_axis, profile.ψ, profile.z; observation_line..., label = "RAPID")
        end
        if !isnothing(index)
            lines!(index_axis, six_month_means(index.year, index.ψ)...; observation_line..., label = "RAPID")
        end
        for (n, output) in enumerate(outputs)
            ψ = atlantic_overturning_streamfunction(output)
            run_lines!(profile_axis, ψ.ψ[rapid_row(output), :], ψ.depths, n, output.label)
            series = rapid_index_series(output)
            dates = record_dates(output, series.years .* seconds_per_year)
            calendar_years = forcing_years(output, series.years .* seconds_per_year) .+ (Dates.dayofyear.(dates) .- 0.5) ./ 365.25
            last_cycle = series.years .≥ maximum(series.years) - (output.forcing_years[2] - output.forcing_years[1] + 1)
            run_lines!(index_axis, six_month_means(calendar_years[last_cycle], series.ψ[last_cycle])..., n, output.label)
        end
        ylims!(profile_axis, -maximum_depth, 0)
    end
end

north_atlantic_map!(position, data; kwargs...) =
    map_panel!(position, map_longitudes, map_latitudes, data; dest = "+proj=stere +lat_0=55 +lon_0=-35 +lat_ts=55",
               longitude_limits = (-80, 10), latitude_limits = (35, 70), kwargs...)

ClimaOcean.barotropic_streamfunction(outputs::OMIPOutput...) =
    maps_figure(outputs, [map_row("barotropic streamfunction", barotropic_streamfunction_map, :balance, (-60, 60), "Sv"),
                          map_row("North Atlantic", barotropic_streamfunction_map, :balance, (-40, 40), "Sv";
                                  draw = north_atlantic_map!)])

const strait_observations = (bering = (mean = 1.0, low = 0.8, high = 1.2),
                             indonesian_throughflow = (mean = -15.0, low = -18.7, high = -10.7))

# Drake Passage transport estimates [Sv]: Chidichimo et al. (2014), Lumpkin and Speer (2007), Whitworth and
# Peterson (1985), Xu et al. (2020) and Koenig et al. (2014)
const drake_observations = [(127.7, 8.1), (129.0, 6.0), (134.0, 13.0), (135.3, 10.2), (141.0, 2.7)]

function ClimaOcean.strait_transports(outputs::OMIPOutput...)
    titles = ("Bering Strait", "Drake Passage", "Indonesian Throughflow")
    names = (:bering, :drake, :indonesian_throughflow)
    panels = [position -> time_axis(position; title, ylabel = "transport (Sv)") for title in titles]

    return lines_figure(panels) do axes
        for (axis, name) in zip(axes, names)
            if haskey(strait_observations, name)
                observed = strait_observations[name]
                hspan!(axis, observed.low, observed.high; observation_band...)
                hlines!(axis, observed.mean; observation_line..., label = "observations")
            end
            for (n, output) in enumerate(outputs)
                series = strait_transport_series(output)
                run_lines!(axis, annual_means(series.years, getproperty(series, name))..., n, output.label)
            end
        end
        last_year = maximum(o -> maximum(strait_transport_series(o).years), outputs)
        x = last_year .* (1.03 .+ 0.02 .* (0:length(drake_observations) - 1))
        errorbars!(axes[2], x, first.(drake_observations), last.(drake_observations); color = observation_color)
        scatter!(axes[2], x, first.(drake_observations); color = observation_color, label = "observations")
    end
end

# Southward export out of the Arctic is negative: de Steur et al. (2009) and Spreen et al. (2009) for Fram,
# Curry et al. (2014) for Davis
const arctic_freshwater_observations = (fram = (liquid = (-2660.0, 500.0), solid = (-2300.0, 500.0)),
                                        davis = (liquid = (-2930.0, 400.0), solid = (-315.0, 150.0)))

function ClimaOcean.arctic_freshwater_transports(outputs::OMIPOutput...)
    return omip_figure(3, 2; column_width = 520, row_height = 400) do figure
        axes = Dict()
        for (r, strait) in enumerate((:fram, :davis)), (c, phase) in enumerate((:liquid, :solid))
            axis = time_axis(figure[r, c]; title = string(uppercasefirst(string(strait)), " Strait, ", phase),
                             ylabel = "freshwater transport (km³ yr⁻¹)")
            center, spread = getproperty(getproperty(arctic_freshwater_observations, strait), phase)
            hspan!(axis, center - spread, center + spread; observation_band...)
            hlines!(axis, center; observation_line..., label = "observations")
            hlines!(axis, 0; color = (:black, 0.3), linewidth = 1)
            axes[(strait, phase)] = axis
        end
        for (n, output) in enumerate(outputs)
            series = arctic_freshwater_series(output)
            for strait in (:fram, :davis)
                run_lines!(axes[(strait, :liquid)], annual_means(series.liquid_years, getproperty(series.liquid, strait))..., n, output.label)
                run_lines!(axes[(strait, :solid)], annual_means(series.solid_years, getproperty(series.solid, strait))..., n, output.label)
            end
        end
        legend!(figure[1:2, 3], axes[(:fram, :liquid)])
    end
end

#####
##### Every figure
#####

const omip_figure_functions = (ClimaOcean.sea_surface_temperature_bias, ClimaOcean.sea_surface_salinity_bias,
                               ClimaOcean.sea_surface_height, ClimaOcean.mixed_layer_depth,
                               ClimaOcean.sea_ice_concentration, ClimaOcean.sea_ice_concentration_bias,
                               ClimaOcean.surface_fluxes, ClimaOcean.wind_stress,
                               ClimaOcean.sea_surface_height_variability, ClimaOcean.sea_ice_extent,
                               ClimaOcean.sea_ice_area, ClimaOcean.arctic_sea_ice_volume,
                               ClimaOcean.sea_ice_area_time_series, ClimaOcean.arctic_sea_ice_volume_time_series,
                               ClimaOcean.kinetic_energy, ClimaOcean.global_drift, ClimaOcean.mean_profiles,
                               ClimaOcean.zonal_mean_sections, ClimaOcean.zonal_mean_drift,
                               ClimaOcean.zonal_mean_mixed_layer_depth, ClimaOcean.drift_profiles,
                               ClimaOcean.strait_transports, ClimaOcean.atlantic_overturning,
                               ClimaOcean.near_surface_currents, ClimaOcean.equatorial_undercurrent,
                               ClimaOcean.atlantic_overturning_at_rapid, ClimaOcean.sea_level_drift,
                               ClimaOcean.content_conservation, ClimaOcean.barotropic_streamfunction,
                               ClimaOcean.arctic_freshwater_transports, ClimaOcean.zonal_rms_bias,
                               ClimaOcean.rms_bias_profiles)

function ClimaOcean.omip_figures(outputs::OMIPOutput...; directory = ".")
    mkpath(directory)
    for (n, figure_function) in enumerate(omip_figure_functions)
        name = string(lpad(n, 2, '0'), '_', nameof(figure_function), ".png")
        @info "Drawing $name"
        save(joinpath(directory, name), figure_function(outputs...))
    end
    return directory
end
