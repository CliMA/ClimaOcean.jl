# Metadata

The [`ClimaOcean.DataWrangling`](@ref ClimaOcean.DataWrangling) module represents each external
dataset with a lightweight metadata object. A piece of metadata does not hold the data itself;
instead, it records everything needed to locate, download, and interpret a variable when the model
asks for it.
For example, to build a "Metadatum" (a single snapshot in time) representing temperature
on January 1st, 2010 in the [EN4 dataset](https://www.metoffice.gov.uk/hadobs/en4/), we write

```@example metadata
using ClimaOcean, Dates

metadatum = Metadatum(:temperature;
                      dataset = EN4Monthly(),
                      date = Date(2010, 1, 1)) 
```

To download and instantiate the data, we use `set!`,

```@example metadata
using Oceananigans

grid = LatitudeLongitudeGrid(size=(360, 90, 1), latitude=(-90, 90), longitude=(0, 360), z=(0, 1))
T = CenterField(grid)
set!(T, metadatum)
```

and then we can plot it:

```@example metadata 
using CairoMakie
heatmap(T)
```

The key ingredients stored in a [`Metadata`](@ref) or [`Metadatum`](@ref) object are

- the variable `Symbol` (for example `:temperature` or `:u_velocity`);
- the concrete dataset type (such as [`EN4Monthly`](@ref), [`ECCO2Daily`](@ref), or [`GLORYSMonthly`](@ref));
- the temporal coverage: either a single timestamp (`Metadatum`) or a range/vector of dates (`Metadata`);
- an optional [`BoundingBox`](@ref) describing regional subsets in longitude, latitude, or depth;
- the on-disk directory where the dataset should be cached.

This bookkeeping lets downstream utilities (for example `set!` or `FieldTimeSeries`)
request exactly the slices of data they need, and it keeps track of where those slices live so we do
not redownload them unnecessarily.

## Supported datasets

ClimaOcean currently ships connectors for the following data products:

- [`ETOPO2022`](@ref) | [Supported variables](@ref dataset-etopo2022-vars) | [NOAA ETOPO 2022 overview](https://www.ncei.noaa.gov/products/etopo-global-relief-model)
- [`ECCO2Monthly`](@ref) | [Supported variables](@ref dataset-ecco2monthly-vars) | [ECCO2 documentation](https://ecco.jpl.nasa.gov/products/all/)
- [`ECCO2Daily`](@ref) | [Supported variables](@ref dataset-ecco2daily-vars) | [ECCO2 documentation](https://ecco.jpl.nasa.gov/products/all/)
- [`ECCO4Monthly`](@ref) | [Supported variables](@ref dataset-ecco4monthly-vars) | [ECCO V4r4 product guide](https://ecco-group.org/products-ECCO-V4r4.htm)
- [`EN4Monthly`](@ref) | [Supported variables](@ref dataset-en4monthly-vars) | [Met Office EN4 overview](https://www.metoffice.gov.uk/hadobs/en4/)
- [`GLORYSDaily`](@ref) | [Supported variables](@ref dataset-glorysdaily-vars) | [Copernicus GLORYS product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description)
- [`GLORYSMonthly`](@ref) | [Supported variables](@ref dataset-glorysmonthly-vars) | [Copernicus GLORYS product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description)
- [`RepeatYearJRA55`](@ref) | [Supported variables](@ref dataset-repeatyearjra55-vars) | [JRA-55 Reanalysis](https://jra.kishou.go.jp/JRA-55/index_en.html)
- [`MultiYearJRA55`](@ref) | [Supported variables](@ref dataset-multiyearjra55-vars) | [JRA-55 Reanalysis](https://jra.kishou.go.jp/JRA-55/index_en.html)

(@id dataset-etopo2022-vars)
### Supported variables for ETOPO2022

- `:bottom_height` - Global bathymetry/topography on a 1 arc-minute grid.

(@id dataset-ecco2monthly-vars)
### Supported variables for ECCO2Monthly

- `:temperature` - Potential temperature on the ECCO2 three-dimensional grid (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:free_surface` - Sea surface height anomaly (m).
- `:sea_ice_thickness` - Effective sea-ice thickness (m).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:net_heat_flux` - Net surface heat flux into the ocean (W m^-2).

(@id dataset-ecco2daily-vars)
### Supported variables for ECCO2Daily

- `:temperature` - Potential temperature on the ECCO2 three-dimensional grid (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:free_surface` - Sea surface height anomaly (m).
- `:sea_ice_thickness` - Effective sea-ice thickness (m).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:net_heat_flux` - Net surface heat flux into the ocean (W m^-2).

(@id dataset-ecco4monthly-vars)
### Supported variables for ECCO4Monthly

- `:temperature` - Potential temperature on the ECCO V4 grid (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:free_surface` - Sea surface height anomaly (m).
- `:sea_ice_thickness` - Effective sea-ice thickness (m).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:net_heat_flux` - Net surface heat flux into the ocean (W m^-2).
- `:sensible_heat_flux` - Surface sensible heat flux (W m^-2).
- `:latent_heat_flux` - Surface latent heat flux (W m^-2).
- `:net_longwave` - Net longwave radiation at the surface (W m^-2).
- `:downwelling_shortwave` - Downward shortwave radiation at the surface (W m^-2).
- `:downwelling_longwave` - Downward longwave radiation at the surface (W m^-2).

(@id dataset-en4monthly-vars)
### Supported variables for EN4Monthly

- `:temperature` - Objective analyses of ocean temperature (degC).
- `:salinity` - Objective analyses of ocean salinity (psu).

(@id dataset-glorysdaily-vars)
### Supported variables for GLORYSDaily

- `:temperature` - Potential temperature (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:sea_ice_thickness` - Sea-ice thickness (m).
- `:sea_ice_u_velocity` - Zonal sea-ice drift (m s^-1).
- `:sea_ice_v_velocity` - Meridional sea-ice drift (m s^-1).
- `:free_surface` - Sea surface height (m).
- `:depth` - Static bathymetry/depth (m).

(@id dataset-glorysmonthly-vars)
### Supported variables for GLORYSMonthly

- `:temperature` - Potential temperature (degC).
- `:salinity` - Practical salinity (psu).
- `:u_velocity` - Zonal velocity component (m s^-1).
- `:v_velocity` - Meridional velocity component (m s^-1).
- `:sea_ice_concentration` - Sea-ice area fraction (dimensionless).
- `:sea_ice_thickness` - Sea-ice thickness (m).
- `:sea_ice_u_velocity` - Zonal sea-ice drift (m s^-1).
- `:sea_ice_v_velocity` - Meridional sea-ice drift (m s^-1).
- `:free_surface` - Sea surface height (m).
- `:depth` - Static bathymetry/depth (m).

(@id dataset-repeatyearjra55-vars)
### Supported variables for RepeatYearJRA55

- `:temperature` - 2 m air temperature (K).
- `:specific_humidity` - 2 m specific humidity (kg kg^-1).
- `:eastward_velocity` - 10 m eastward wind (m s^-1).
- `:northward_velocity` - 10 m northward wind (m s^-1).
- `:sea_level_pressure` - Sea-level pressure (Pa).
- `:downwelling_shortwave_radiation` - Downward shortwave radiation (W m^-2).
- `:downwelling_longwave_radiation` - Downward longwave radiation (W m^-2).
- `:rain_freshwater_flux` - Liquid precipitation flux (kg m^-2 s^-1).
- `:snow_freshwater_flux` - Solid precipitation flux (kg m^-2 s^-1).
- `:river_freshwater_flux` - River discharge flux (kg m^-2 s^-1).
- `:iceberg_freshwater_flux` - Iceberg calving flux (kg m^-2 s^-1).

(@id dataset-multiyearjra55-vars)
### Supported variables for MultiYearJRA55

- `:temperature` - 2 m air temperature (K).
- `:specific_humidity` - 2 m specific humidity (kg kg^-1).
- `:eastward_velocity` - 10 m eastward wind (m s^-1).
- `:northward_velocity` - 10 m northward wind (m s^-1).
- `:sea_level_pressure` - Sea-level pressure (Pa).
- `:downwelling_shortwave_radiation` - Downward shortwave radiation (W m^-2).
- `:downwelling_longwave_radiation` - Downward longwave radiation (W m^-2).
- `:rain_freshwater_flux` - Liquid precipitation flux (kg m^-2 s^-1).
- `:snow_freshwater_flux` - Solid precipitation flux (kg m^-2 s^-1).
- `:river_freshwater_flux` - River discharge flux (kg m^-2 s^-1).
- `:iceberg_freshwater_flux` - Iceberg calving flux (kg m^-2 s^-1).


### Metadata helpers used below

```@setup metadata
using ClimaOcean
using ClimaOcean.DataWrangling: retrieve_data, download_dataset, metadata_path
using Dates
using Oceananigans
using Oceananigans.Units: day
using CairoMakie
using CairoMakie: Colorbar, Observable, record

const metadata_assets_dir = joinpath(@__DIR__, "assets")
isdir(metadata_assets_dir) || mkpath(metadata_assets_dir)

finite_extrema(A) = begin
    values = filter(!isnan, vec(Array(A)))
    isempty(values) ? (0.0, 1.0) : extrema(values)
end
```

The setup above loads the functions used throughout this tutorial, defines a small helper for
computing color limits that ignore `NaN`s, and ensures we have a place to save any animations.

## Single-date snapshots

[`Metadatum`](@ref) is a convenience constructor for a `Metadata` object that points to a single
timestamp. It is ideal for one-off snapshots, for example to initialize a field or inspect a
particular profile.


The metadata above is lazy—it simply remembers that we want EN4 temperature for January 2010. When
we pass it to `set!` the first time, ClimaOcean downloads (or reuses) the corresponding NetCDF file,
interpolates onto the requested grid, and fills the target field.

```@example metadata
Nx, Ny, Nz = 360, 180, 4
underlying_grid = LatitudeLongitudeGrid(size = (Nx, Ny, Nz),
                                        longitude = (0, 360),
                                        latitude = (-90, 90),
                                        z = (-5000, 0))

bottom_height = regrid_bathymetry(underlying_grid)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

T = CenterField(grid)
set!(T, snapshot)  # downloads on first use and fills the field

parent(T) |> size
```

The field now contains the requested snapshot. Plotting several depth levels is a convenient way to
check the structure of the imported data.

```@example metadata
levels = (grid.Nz, max(grid.Nz - 1, 1), 1)  # bottom, near-bottom, and surface
temperature_slices = [Array(view(T, :, :, level)) for level in levels]
color_limits = finite_extrema(parent(T))

fig = Figure(resolution = (960, 320))
heatmaps = Vector{Any}()

for (col, (level, slice)) in enumerate(zip(levels, temperature_slices))
    ax = Axis(fig[1, col],
              title = "Level $level",
              xlabel = col == 1 ? "i-index" : "",
              ylabel = col == 1 ? "j-index" : "")
    hm = heatmap!(ax, slice;
                  colormap = :thermal,
                  colorrange = color_limits,
                  nan_color = :lightgray)
    push!(heatmaps, hm)
end

Colorbar(fig[1, length(levels) + 1],
         heatmaps[1],
         label = "Temperature (°C)")

fig
```

If you prefer to flip through every vertical level, the short snippet below records a movie. The
resulting `mp4` is saved in `docs/src/assets/` so it can be embedded or inspected independently.

```@example metadata
all_levels = collect(grid.Nz:-1:1)

slice_at(level) = Array(view(T, :, :, level))
temperature_data = Observable(slice_at(first(all_levels)))
color_limits = finite_extrema(parent(T))

fig = Figure(resolution = (520, 360))
ax = Axis(fig[1, 1],
          xlabel = "i-index",
          ylabel = "j-index",
          title = "Depth level $(first(all_levels))")
hm = heatmap!(ax, temperature_data;
              colormap = :thermal,
              colorrange = color_limits,
              nan_color = :lightgray)
Colorbar(fig[1, 2], hm, label = "Temperature (°C)")

animation_path = joinpath(metadata_assets_dir, "en4_temperature_levels.mp4")
record(fig, animation_path, all_levels; framerate = 2) do level
    temperature_data[] = slice_at(level)
    ax.title = "Depth level $level"
end

println("Saved animation: $(animation_path)")
fig
```

## ECCO time series and derived diagnostics

`Metadata` objects become especially useful when we want a time series. Below we grab daily surface
velocities from ECCO, compute a quick-and-dirty relative vorticity diagnostic, and save both a
figure and an animation of the evolving field. This example assumes that you have configured the
ECCO WebDAV credentials described in the [`ECCO README`](https://github.com/CliMA/ClimaOcean.jl/blob/main/src/DataWrangling/ECCO/README.md).

!!! warning "Credentials required"
    Set the `ECCO_USERNAME` and `ECCO_WEBDAV_PASSWORD` environment variables before running the
    example. The code below skips execution if the credentials are missing.

```@example metadata
has_ecco_credentials =
    !isnothing(get(ENV, "ECCO_USERNAME", nothing)) &&
    !isnothing(get(ENV, "ECCO_WEBDAV_PASSWORD", nothing))

if has_ecco_credentials
    start_date = DateTime(1994, 1, 1)
    end_date   = DateTime(1994, 1, 5)
    region     = BoundingBox(longitude = (220, 250), latitude = (-10, 10))

    u_series = Metadata(:u_velocity;
                        dataset = ECCO2Daily(),
                        start_date,
                        end_date,
                        bounding_box = region)
    v_series = Metadata(:v_velocity;
                        dataset = ECCO2Daily(),
                        start_date,
                        end_date,
                        bounding_box = region)

    download_dataset(u_series)
    download_dataset(v_series)

    using NCDatasets

    sample_path = metadata_path(first(u_series))
    ds = Dataset(sample_path)
    lon_full = ds["XC"][:, 1]
    lat_full = ds["YC"][1, :]
    close(ds)

    lon_indices = findall(λ -> region.longitude[1] <= λ <= region.longitude[2], lon_full)
    lat_indices = findall(φ -> region.latitude[1]  <= φ <= region.latitude[2], lat_full)
    lon_axis = lon_full[lon_indices]
    lat_axis = lat_full[lat_indices]

    earth_radius = 6.371e6

    function surface_vorticity(u_slice, v_slice, lon_axis, lat_axis)
        Nx, Ny = size(u_slice)
        ζ = zeros(Float32, Nx, Ny)
        lon_rad = deg2rad.(lon_axis)
        lat_rad = deg2rad.(lat_axis)
        cos_lat = cos.(lat_rad)

        for j in 2:Ny-1
            Δy = earth_radius * (lat_rad[j+1] - lat_rad[j-1]) / 2
            for i in 2:Nx-1
                Δx = earth_radius * cos_lat[j] * (lon_rad[i+1] - lon_rad[i-1]) / 2
                ζ[i, j] = (v_slice[i+1, j] - v_slice[i-1, j]) / (2Δx) -
                          (u_slice[i, j+1] - u_slice[i, j-1]) / (2Δy)
            end
            ζ[1, j]  = ζ[2, j]
            ζ[end, j] = ζ[end-1, j]
        end

        ζ[:, 1]  .= ζ[:, 2]
        ζ[:, end] .= ζ[:, end-1]
        return ζ
    end

    vorticity_frames = Matrix{Float32}[]
    frame_dates = DateTime[]

    for (u_meta, v_meta) in zip(u_series, v_series)
        u_data = retrieve_data(u_meta)
        v_data = retrieve_data(v_meta)

        surface_u = Array(u_data[lon_indices, lat_indices, end])
        surface_v = Array(v_data[lon_indices, lat_indices, end])

        push!(vorticity_frames, surface_vorticity(surface_u, surface_v, lon_axis, lat_axis))
        push!(frame_dates, u_meta.dates)
    end

    values = Float32[]
    for frame in vorticity_frames
        append!(values, filter(!isnan, vec(frame)))
    end
    color_limits = isempty(values) ? (-1f-6, 1f-6) : extrema(values)

    vorticity_data = Observable(first(vorticity_frames))
    fig = Figure(resolution = (720, 400))
    ax = Axis(fig[1, 1],
              xlabel = "Longitude (°E)",
              ylabel = "Latitude (°N)",
              title = "Surface vorticity on $(first(frame_dates))")
    hm = heatmap!(ax, lon_axis, lat_axis, vorticity_data;
                  colormap = :balance,
                  colorrange = color_limits,
                  nan_color = :lightgray)
    Colorbar(fig[1, 2], hm, label = "s⁻¹")

    animation_path = joinpath(metadata_assets_dir, "ecco_surface_vorticity.mp4")
    record(fig, animation_path, eachindex(vorticity_frames); framerate = 2) do idx
        vorticity_data[] = vorticity_frames[idx]
        ax.title = "Surface vorticity on $(frame_dates[idx])"
    end

    println("Saved animation: $(animation_path)")
    fig
else
    println("Skipping ECCO example because ECCO credentials are not configured.")
end
```

The animation illustrates how convenient it is to iterate over `Metadata`—each iteration produces a
[`Metadatum`](@ref) for the next time slice, so film loops naturally build out of plain Julia code.

## CopernicusMarine subsets with bounding boxes

CopernicusMarine's GLORYS products expose high-resolution data, so bounding boxes are essential
when exploring a sub-region. When the optional
[`CopernicusMarine.jl`](https://github.com/CliMA/CopernicusMarine.jl) extension is available,
`download_dataset` automatically asks the Copernicus subset API for the specified spatial window.

!!! warning "Copernicus login required"
    Export `COPERNICUS_USERNAME` and `COPERNICUS_PASSWORD` to allow ClimaOcean to query the Copernicus
    Marine Service. The example below skips execution if the credentials are missing.

```@example metadata
has_cmems_credentials =
    !isnothing(get(ENV, "COPERNICUS_USERNAME", nothing)) &&
    !isnothing(get(ENV, "COPERNICUS_PASSWORD", nothing))

if has_cmems_credentials
    bbox = BoundingBox(longitude = (210, 250),
                       latitude = (38, 52),
                       z = (-2000, -50))

    copernicus_snapshot = Metadatum(:temperature;
                                    dataset = GLORYSMonthly(),
                                    date = DateTime(2018, 7, 1),
                                    bounding_box = bbox)

    download_dataset(copernicus_snapshot)
    native_field = Field(copernicus_snapshot)
    native_grid = native_field.grid

    function axis_from_grid(grid)
        lon = [Oceananigans.Grids.node(i, 1, 1, grid, Center(), Center(), Center())[1]
               for i in 1:grid.Nx]
        lat = [Oceananigans.Grids.node(1, j, 1, grid, Center(), Center(), Center())[2]
               for j in 1:grid.Ny]
        return lon, lat
    end

    lon_axis, lat_axis = axis_from_grid(native_grid)
    surface_slice = Array(view(native_field, :, :, native_grid.Nz))
    color_limits = finite_extrema(surface_slice)

    fig = Figure(resolution = (640, 420))
    ax = Axis(fig[1, 1],
              xlabel = "Longitude (°E)",
              ylabel = "Latitude (°N)",
              title = "GLORYS temperature on $(copernicus_snapshot.dates)")
    hm = heatmap!(ax, lon_axis, lat_axis, surface_slice;
                  colormap = :thermal,
                  colorrange = color_limits,
                  nan_color = :lightgray)
    Colorbar(fig[1, 2], hm, label = "Temperature (°C)")

    println("Native grid size: $(size(surface_slice)) (lon × lat)")
    fig
else
    println("Skipping Copernicus example because COPERNICUS credentials are not configured.")
end
```

Notice how the grid dimensions reflect the geographic subset; ClimaOcean trims the native grid to
match the longitude, latitude, and depth bounds specified by the bounding box. This keeps downloads
lightweight and makes it easy to work on regional problems.

## Dataset-based restoring as a forcing

Metadata does more than populate initial conditions. The [`DatasetRestoring`](@ref) forcing wraps a
`Metadata` time series, handles downloads and caching, and continuously nudges a simulation toward
the prescribed data. You can interpolate on the fly or pre-interpolate on a user-supplied grid,
choose inpainting strategies, and mask regions where restoring should be weak or absent.

The snippet below builds an EN4 temperature restoring term that relaxes a model toward three months
of climatology on a 72 × 36 × 10 grid, tapering the tendency toward the poles.

```@example metadata
restoring_metadata = Metadata(:temperature;
                              dataset = EN4Monthly(),
                              start_date = Date(2010, 1, 1),
                              end_date   = Date(2010, 3, 1))

restoring_grid = LatitudeLongitudeGrid(size = (72, 36, 10),
                                       longitude = (0, 360),
                                       latitude = (-80, 80),
                                       z = (-2000, 0),
                                       halo = (3, 3, 3))

polar_mask = LinearlyTaperedPolarMask(northern = (70, 75),
                                      southern = (-75, -70),
                                      z = (-100, 0))

temperature_restoring = DatasetRestoring(restoring_metadata,
                                         restoring_grid;
                                         rate = 1 / (10day),
                                         mask = polar_mask)

fields = (T = CenterField(restoring_grid),)
fill!(fields.T, 0)
clock = Clock(; time = 0.0)
sample_tendency = temperature_restoring(10, 10, 5,
                                        restoring_grid,
                                        clock,
                                        fields)
println("Sample restoring tendency (s⁻¹ × native anomaly): $(sample_tendency)")

temperature_restoring
```

Attach `temperature_restoring` to an `Oceananigans` model just like any other forcing—for example,
`model.forcing.T = temperature_restoring`. The helper writes the necessary `FieldTimeSeries` to
cache, manages on-disk horizontal masks produced by inpainting, and keeps the lookup asynchronous
so the simulation proceeds smoothly.
