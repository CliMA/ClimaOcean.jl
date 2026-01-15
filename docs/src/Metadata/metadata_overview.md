# Metadata

[`Metadata`](@ref) is an abstraction that _represents_ data, but does not embody it.
Unlike [`Oceananigans.Field`](https://clima.github.io/OceananigansDocumentation/stable/appendix/library/#Oceananigans.Fields.Field-Union{Tuple{Oceananigans.Grids.AbstractGrid},%20Tuple{LZ},%20Tuple{LY},%20Tuple{LX},%20Tuple{Oceananigans.Grids.AbstractGrid,%20DataType}}%20where%20{LX,%20LY,%20LZ}),
which points to an array occupying space in memory,
`Metadata` only contains information about where files are stored, their origin, the grid
they live on, and the date(s) they correspond to (if any).

```@docs
Metadata
```

When `Metadata` represents just one date, we call it [`Metadatum`](@ref).
For example, consider global temperature from January 1st, 2010 from the
[EN4 dataset](https://www.metoffice.gov.uk/hadobs/en4/),

```@example metadata
using ClimaOcean, Dates

metadatum = Metadatum(:temperature;
                      dataset = EN4Monthly(),
                      date = Date(2010, 1, 1))
```

To materialize the data described by a `metadatum`, we wrap it in an Oceananigans' `Field`,

```@example metadata
using Oceananigans

T_native = Field(metadatum)
```

We can also interpolate the data on a user-defined grid by using the function `set!`,

```@example metadata
grid = LatitudeLongitudeGrid(size = (360, 90, 1),
                             latitude = (-90, 90),
                             longitude = (0, 360),
                             z = (0, 1))
T = CenterField(grid)
set!(T, metadatum)
```

and then we can plot it:

```@example metadata
using CairoMakie
heatmap(T)
```

This looks a bit odd, but less so if we download bathymetry (for which we also use `Metadata`
under the hood) to create a temperature field with a land mask,

```@example metadata
bottom_height = regrid_bathymetry(grid)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))
T = CenterField(grid)
set!(T, metadatum)
heatmap(T)
```

The key ingredients stored in a [`Metadata`](@ref) or [`Metadatum`](@ref) object are:

- the variable name (for example `:temperature` or `:u_velocity`);
- the dataset (such as `EN4Monthly`, `ECCO2Daily`, or `GLORYSMonthly`);
- the temporal coverage: either a single timestamp (`Metadatum`) or a range/vector of dates (`Metadata`);
- an optional [`BoundingBox`](@ref ClimaOcean.DataWrangling.BoundingBox) describing regional subsets in
  longitude, latitude, or depth;
- the on-disk `dir`ectory where the dataset are be cached.

This bookkeeping lets downstream utilities (for example `set!` or `FieldTimeSeries`) request exactly the
slices of data they need, and it keeps track of where those slices live so we do not redownload
them unnecessarily.

## Supported datasets

ClimaOcean currently ships connectors for the following data products:

| Dataset            | Supported Variables                                      | Documentation Link                                                                                 |
|--------------------|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------------|
| `ETOPO2022`        | [Supported variables](@ref dataset-etopo2022-vars)        | [NOAA ETOPO 2022 overview](https://www.ncei.noaa.gov/products/etopo-global-relief-model)           |
| `ECCO2Monthly`     | [Supported variables](@ref dataset-ecco2monthly-vars)     | [ECCO2 documentation](https://ecco.jpl.nasa.gov/products/all/)                                     |
| `ECCO2Daily`       | [Supported variables](@ref dataset-ecco2daily-vars)       | [ECCO2 documentation](https://ecco.jpl.nasa.gov/products/all/)                                     |
| `ECCO4Monthly`     | [Supported variables](@ref dataset-ecco4monthly-vars)     | [ECCO V4r4 product guide](https://ecco-group.org/products-ECCO-V4r4.htm)                           |
| `EN4Monthly`       | [Supported variables](@ref dataset-en4monthly-vars)       | [Met Office EN4 overview](https://www.metoffice.gov.uk/hadobs/en4/)                                |
| `GLORYSDaily`      | [Supported variables](@ref dataset-glorysdaily-vars)      | [Copernicus GLORYS product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) |
| `GLORYSMonthly`    | [Supported variables](@ref dataset-glorysmonthly-vars)    | [Copernicus GLORYS product page](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) |
| `RepeatYearJRA55`  | [Supported variables](@ref dataset-repeatyearjra55-vars)  | [JRA-55 Reanalysis](https://www.data.jma.go.jp/jra/html/JRA-55/index_en.html)                                 |
| `MultiYearJRA55`   | [Supported variables](@ref dataset-multiyearjra55-vars)   | [JRA-55 Reanalysis](https://www.data.jma.go.jp/jra/html/JRA-55/index_en.html)                                 |
