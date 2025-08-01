# Datasets

The [`DataWrangling`](@ref ClimaOcean.DataWrangling) module allows us to use various datasets
in ClimaOcean simulations.

A central concept in the `DataWrangling` module is that of [`Metadata`](@ref).
Metadata is a collection of information for the dataset, like the dataset's name, the variable
that it involves, the date range, where the dataset is stored locally.

For example, below we construct a metadata for the temperature variable from the EN4Monthly dataset.

```@example metadata
using ClimaOcean
using Dates

T_meta = Metadata(:temperature, dataset=EN4Monthly())
```

We can use `start_date` and `end_date` keyword arguments to select only some of the data
available in the dataset. For example, to get only the first quarter of 2010 we do

```@example metadata
T_meta = Metadata(:temperature, dataset=EN4Monthly(),
                  start_date = Date(2010, 1),
                  end_date = Date(2010, 3))
```

We can use `set!` on a field or a model

```@example metadata
using Oceananigans

Nx, Ny, Nz = 720, 360, 4
underlying_grid = LatitudeLongitudeGrid(size=(Nx, Ny, Nz),
                                        longitude = (0, 360),
                                        latitude = (-90, 90),
                                        z = (-5000, 0))

bottom_height = regrid_bathymetry(underlying_grid)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

T = CenterField(grid)
```

Note that `T` field above is empty. Now we can set it to the first time slice of `T_meta`.

```@example metadata
set!(T, first(T_meta))
```

Now `T` has values. Also note how `set!` above triggered downloading the dataset that
corresponds to `T_meta`.

```@example metadata
using CairoMakie

heatmap(view(T, :, :, grid.Nz))
```
