#####
##### Reading the output of a run, once
#####

cached(f, output::OMIPOutput, key) = get!(f, output.cache, key)

timeseries(output, group, name) = cached(output, Symbol(group, '/', name)) do
    FieldTimeSeries(output_path(output, group), string(name); backend = OnDisk())
end

function first_output_file(output, group)
    path = output_path(output, group)
    isfile(path) && return path
    stem = basename(path)[1:end-5]
    return joinpath(output.directory, stem * "_part1.jld2")
end

has_timeseries(output, group, name) =
    jldopen(file -> haskey(file, "timeseries/$name"), first_output_file(output, group))

model_grid(output) = cached(() -> timeseries(output, "surface", "tos").grid, output, :grid)

const seconds_per_year = 365.25 * 86400

in_years(times) = times ./ seconds_per_year

function window_records(output, times)
    t₁, t₂ = output.averaging_window
    records = findall(t -> t₁ ≤ t ≤ t₂, times)
    isempty(records) && error("$(output.label) has no record in the averaging window $(output.averaging_window)")
    return records
end

record_array(fts, n) = Array(interior(fts[n]))

"""
The mean of `name` in `group` over the records of the averaging window.
"""
time_mean(output, group, name) = cached(output, Symbol(:time_mean, '/', group, '/', name)) do
    fts = timeseries(output, group, name)
    records = window_records(output, fts.times)
    total = record_array(fts, first(records))
    for n in records[2:end]
        total .+= record_array(fts, n)
    end
    total ./ length(records)
end

# The calendar date at the middle of each averaging interval
function record_dates(output, times)
    half_interval = length(times) > 1 ? (times[2] - times[1]) / 2 : 0.0
    return [output.start_date + Millisecond(round(Int, 1000 * (t - half_interval))) for t in times]
end

"""
The mean of `name` in `group` over the records of each calendar month of the averaging window.
"""
monthly_means(output, group, name) = cached(output, Symbol(:monthly_means, '/', group, '/', name)) do
    fts = timeseries(output, group, name)
    records = window_records(output, fts.times)
    months = month.(record_dates(output, fts.times[records]))
    map(1:12) do m
        in_month = records[months .== m]
        isempty(in_month) && error("$(output.label) has no record in month $m of its averaging window")
        sum(record_array(fts, n) for n in in_month) ./ length(in_month)
    end
end

# The forcing years a run went through, as calendar years of the observations
function forcing_years(output, times)
    first_year, last_year = output.forcing_years
    cycle = last_year - first_year + 1
    simulated_years = year.(record_dates(output, times)) .- year(output.start_date)
    return first_year .+ mod.(simulated_years, cycle)
end

function averaged_forcing_years(output)
    times = timeseries(output, "surface", "siconc").times
    return extrema(forcing_years(output, times[window_records(output, times)]))
end

#####
##### Masks and coordinates of the model grid
#####

ocean_mask(output) = cached(output, :ocean_mask) do
    grid = model_grid(output)
    Nx, Ny, Nz = size(grid)
    [!inactive_node(i, j, k, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
end

surface_mask(output) = ocean_mask(output)[:, :, end]

depths(output) = cached(() -> collect(znodes(model_grid(output), Center())), output, :depths)

center_latitudes(output) = cached(output, :center_latitudes) do
    grid = model_grid(output)
    Nx, Ny, _ = size(grid)
    [φnode(i, j, 1, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
end

center_longitudes(output) = cached(output, :center_longitudes) do
    grid = model_grid(output)
    Nx, Ny, _ = size(grid)
    [mod(λnode(i, j, 1, grid, Center(), Center(), Center()), 360) for i in 1:Nx, j in 1:Ny]
end

masked(data, mask) = ifelse.(mask, data, NaN)

"""
The angle between the grid's `i` direction and east at every cell center, as `(cos θ, sin θ)`.
"""
grid_rotation(output) = cached(output, :grid_rotation) do
    grid = model_grid(output)
    Nx, Ny, _ = size(grid)
    cosθ = ones(Nx, Ny)
    sinθ = zeros(Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        Δλ = mod(λnode(i + 1, j, 1, grid, Face(), Center(), Center()) -
                 λnode(i,     j, 1, grid, Face(), Center(), Center()) + 180, 360) - 180
        Δφ = φnode(i + 1, j, 1, grid, Face(), Center(), Center()) - φnode(i, j, 1, grid, Face(), Center(), Center())
        east = cosd(φnode(i, j, 1, grid, Center(), Center(), Center())) * Δλ
        r = hypot(east, Δφ)
        r > 0 && (cosθ[i, j] = east / r; sinθ[i, j] = Δφ / r)
    end
    (cosθ, sinθ)
end

# Grid-aligned components at cell centers rotated to eastward and northward
function eastward_northward(output, u, v)
    cosθ, sinθ = grid_rotation(output)
    return u .* cosθ .- v .* sinθ, u .* sinθ .+ v .* cosθ
end

# Velocity components averaged from their faces to the cell centers, through the fold of the grid
function centered_velocities(output, u_data, v_data)
    grid = model_grid(output)
    Nz = size(u_data, 3)
    LZ = Nz == 1 ? Nothing : Center
    u = Field{Face, Center, LZ}(grid; boundary_conditions = vector_component_boundary_conditions(grid, (Face(), Center(), LZ())))
    v = Field{Center, Face, LZ}(grid; boundary_conditions = vector_component_boundary_conditions(grid, (Center(), Face(), LZ())))
    set!(u, u_data)
    set!(v, v_data)
    fill_halo_regions!(u)
    fill_halo_regions!(v)
    uᶜ = Field(@at((Center, Center, LZ), u))
    vᶜ = Field(@at((Center, Center, LZ), v))
    return Array(interior(compute!(uᶜ))), Array(interior(compute!(vᶜ)))
end
