#####
##### The visual language shared by every figure
#####

# Okabe-Ito extended with Tol (2021): distinguishable under the common forms of colour blindness
const run_palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
                     "#882255", "#117733", "#44AA99", "#332288", "#661100", "#999999"]

run_color(n) = run_palette[mod1(n, length(run_palette))]

const run_linewidth = 2
const observation_color = RGBf(0.45, 0.45, 0.45)
const observation_line = (; color = observation_color, linewidth = 2.5, linestyle = :dash)
const observation_band = (; color = (observation_color, 0.2))

const land_color = RGBf(0.82, 0.82, 0.82)
const month_ticks = (1:12, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"])
const maximum_depth = 5500

const panel_width = 620
const map_height = 360

omip_theme() = Theme(fontsize = 15,
                     Axis = (xgridvisible = false, ygridvisible = false, titlefont = :regular),
                     Legend = (framevisible = false,),
                     Colorbar = (height = 12, ticklabelsize = 13))

function omip_figure(f, columns, rows; column_width = panel_width, row_height = map_height)
    return with_theme(omip_theme()) do
        figure = Figure(size = (column_width * columns + 120, row_height * rows + 80))
        f(figure)
        figure
    end
end

# GeoAxis expects longitudes in [-180, 180)
function centered_longitudes(longitudes, data)
    shifted = mod.(longitudes .+ 180, 360) .- 180
    order = sortperm(shifted)
    return shifted[order], data[order, :]
end

map_levels(colorrange; length = 21) = range(colorrange[1], colorrange[2]; length)

"""
A filled contour map of `data` on `longitudes` and `latitudes` in the projection `dest`, over gray land and
coastlines. Returns the plot so that a row of panels can share one colorbar.
"""
function map_panel!(position, longitudes, latitudes, data;
                    title = "", colormap = :balance, colorrange, dest = "+proj=robin",
                    longitude_limits = nothing, latitude_limits = nothing,
                    edge = nothing, edge_level = 0.15)

    λ, values = centered_longitudes(longitudes, data)
    axis = GeoAxis(position; dest, title, titlegap = 12, xgridvisible = false, ygridvisible = false,
                   xticklabelsvisible = false, yticklabelsvisible = false)

    isnothing(longitude_limits) || xlims!(axis, longitude_limits...)
    isnothing(latitude_limits) || ylims!(axis, latitude_limits...)

    poly!(axis, GeoMakie.land(); color = land_color, strokewidth = 0)
    plot = contourf!(axis, λ, latitudes, values; colormap, levels = map_levels(colorrange),
                     extendlow = :auto, extendhigh = :auto)
    lines!(axis, GeoMakie.coastlines(); color = :black, linewidth = 0.4)

    if !isnothing(edge)
        _, edge_values = centered_longitudes(longitudes, edge)
        filled = map(v -> isfinite(v) ? v : zero(v), edge_values)
        any(<(edge_level), filled) && any(>(edge_level), filled) &&
            contour!(axis, λ, latitudes, filled; levels = [edge_level], color = :black, linewidth = 2)
    end

    return plot
end

global_map!(position, data; kwargs...) = map_panel!(position, map_longitudes, map_latitudes, data; kwargs...)

function polar_map!(position, data; hemisphere, edge = nothing, kwargs...)
    north = hemisphere == :north
    rows = north ? findall(≥(45), map_latitudes) : findall(≤(-45), map_latitudes)
    dest = north ? "+proj=stere +lat_0=90 +lat_ts=70" : "+proj=stere +lat_0=-90 +lat_ts=-70"
    latitude_limits = north ? (45, 90) : (-90, -45)
    edge_rows = isnothing(edge) ? nothing : edge[:, rows]
    return map_panel!(position, map_longitudes, map_latitudes[rows], data[:, rows];
                      dest, latitude_limits, longitude_limits = (-180, 180), edge = edge_rows, kwargs...)
end

# Rows of panels sit at odd rows of the layout, each followed by the colorbar the row shares
panel(figure, row, column) = figure[2row - 1, column]

function row_colorbar!(figure, row, columns, plot, label)
    Colorbar(figure[2row, 1:columns], plot; label, vertical = false, flipaxis = false,
             width = Relative(0.5), tellwidth = false)
    rowgap!(figure.layout, 2row - 1, 4)
    return nothing
end

"""
A latitude--depth or longitude--depth section of `data`, filled between `levels` with thin contour lines on top.
"""
function section_panel!(position, x, z, data; title, levels, colormap = :balance, xlabel = "Latitude",
                        depth_limit = maximum_depth, overlay = nothing, overlay_levels = levels)
    axis = Axis(position; title, xlabel, ylabel = "Depth (m)", backgroundcolor = land_color)
    plot = contourf!(axis, x, z, data; levels, colormap, extendlow = :auto, extendhigh = :auto)
    contour!(axis, x, z, data; levels, color = (:black, 0.5), linewidth = 0.5)
    isnothing(overlay) || contour!(axis, x, z, overlay; levels = overlay_levels, color = :white,
                                   linestyle = :dash, linewidth = 1.2)
    ylims!(axis, -depth_limit, 0)
    return plot
end

"""
A depth--time panel of `data[n, k]` split at 1000 m, so that the upper ocean gets as much room as the abyss.
"""
function depth_time_panel!(position, time_in_years, z, data; title, levels, colormap = :balance)
    layout = GridLayout(position)
    upper = Axis(layout[1, 1]; title, ylabel = "Depth (m)")
    deep = Axis(layout[2, 1]; xlabel = "Year", ylabel = "Depth (m)")
    plot = contourf!(upper, time_in_years, z, data; levels, colormap, extendlow = :auto, extendhigh = :auto)
    contourf!(deep, time_in_years, z, data; levels, colormap, extendlow = :auto, extendhigh = :auto)
    ylims!(upper, -1000, 0)
    ylims!(deep, -maximum_depth, -1000)
    linkxaxes!(upper, deep)
    hidexdecorations!(upper; ticks = false)
    rowgap!(layout, 0)
    return plot
end

time_axis(position; title, ylabel) = Axis(position; title, ylabel, xlabel = "Year")
climatology_axis(position; title, ylabel) = Axis(position; title, ylabel, xlabel = "Month", xticks = month_ticks)

function observation_lines!(axis, x, mean, spread, label)
    band!(axis, x, mean .- spread, mean .+ spread; observation_band...)
    lines!(axis, x, mean; observation_line..., label)
    return nothing
end

run_lines!(axis, x, y, n, label) = lines!(axis, x, y; color = run_color(n), linewidth = run_linewidth, label)

legend!(position, axis) = Legend(position, axis; merge = true, unique = true)
