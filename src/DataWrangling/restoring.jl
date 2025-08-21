using Oceananigans: location
using Oceananigans.Grids: node
using Oceananigans.Fields: interpolate, instantiated_location
using Oceananigans.OutputReaders: Cyclical
using Oceananigans.Utils: Time
using Oceananigans.Architectures: AbstractArchitecture

using JLD2
using NCDatasets

using Dates: Second

import ClimaOcean: stateindex
import Oceananigans.Forcings: regularize_forcing

# Variable names for restorable data
struct Temperature end
struct Salinity end
struct UVelocity end
struct VVelocity end

const oceananigans_fieldnames = Dict(:temperature => Temperature(),
                                     :salinity    => Salinity(),
                                     :u_velocity  => UVelocity(),
                                     :v_velocity  => VVelocity())

@inline Base.getindex(fields, i, j, k, ::Temperature) = @inbounds fields.T[i, j, k]
@inline Base.getindex(fields, i, j, k, ::Salinity)    = @inbounds fields.S[i, j, k]
@inline Base.getindex(fields, i, j, k, ::UVelocity)   = @inbounds fields.u[i, j, k]
@inline Base.getindex(fields, i, j, k, ::VVelocity)   = @inbounds fields.v[i, j, k]

Base.summary(::Temperature) = "temperature"
Base.summary(::Salinity)    = "salinity"
Base.summary(::UVelocity)   = "u_velocity"
Base.summary(::VVelocity)   = "v_velocity"

struct DatasetRestoring{FTS, G, M, V, N}
    field_time_series :: FTS
    native_grid :: G
    mask :: M
    variable_name :: V
    rate :: N
end

Adapt.adapt_structure(to, p::DatasetRestoring) =
    DatasetRestoring(Adapt.adapt(to, p.field_time_series),
                     Adapt.adapt(to, p.native_grid),
                     Adapt.adapt(to, p.mask),
                     Adapt.adapt(to, p.variable_name),
                     Adapt.adapt(to, p.rate))

@inline function (p::DatasetRestoring)(i, j, k, grid, clock, fields)

    # Figure out all the inputs: time, location, and node
    time = Time(clock.time)
    loc = location(p.field_time_series)

    # Possibly interpolate dataset's data from their native grid to simulation grid.
    # Otherwise, simply extract the pre-interpolated data from p.field_time_series.
    if p.native_grid isa Nothing
        ψ_dataset = @inbounds p.field_time_series[i, j, k, time]
    else
        ψ_dataset = interpolate_to_grid(p.field_time_series, i, j, k, p.native_grid, grid, time)
    end

    ψ = @inbounds fields[i, j, k, p.variable_name]
    μ = stateindex(p.mask, i, j, k, grid, clock.time, loc)
    r = p.rate

    return r * μ * (ψ_dataset - ψ)
end

@inline function interpolate_to_grid(fts, i, j, k, native_grid, grid, time)
    times = fts.times
    data = fts.data
    time_indexing = fts.time_indexing
    backend = fts.backend
    loc = instantiated_location(fts)
    X = node(i, j, k, grid, loc...)

    # Interpolate field time series data onto the current node and time
    return interpolate(X, time, data, loc, native_grid, times, backend, time_indexing)
end

default_time_indices_in_memory(metadata) = min(2, length(metadata))
default_time_indices_in_memory(::Metadatum) = 1

"""
    DatasetRestoring(metadata::Metadata,
                     arch_or_grid = CPU();
                     rate,
                     mask = 1,
                     time_indices_in_memory = 2, # Not more than this if we want to use GPU!
                     time_indexing = Cyclical(),
                     inpainting = NearestNeighborInpainting(Inf),
                     cache_inpainted_data = true)

Return a forcing term that restores to data from a dataset. The restoring is applied
as a forcing on the right hand side of the evolution equations, calculated as:

```math
F_ψ = r μ (ψ_{dataset} - ψ)
```

where ``μ`` is the mask, ``r`` is the restoring rate, ``ψ`` is the simulation variable,
and ``ψ_{dataset}`` is the dataset variable that is linearly interpolated in space and time
from the dataset of choice to the simulation grid and time.

Arguments
=========

- `metadata`: The medatada for a dataset variable to restore. Choices for variables include:
  * `:temperature`,
  * `:salinity`,
  * `:u_velocity`,
  * `:v_velocity`,
  * `:sea_ice_thickness`,
  * `:sea_ice_area_fraction`.

- `arch_or_grid`: Either the architecture of the simulation, or a grid on which the data
                  is pre-interpolated when loaded. If an `arch`itecture is provided, such as
                  `arch_or_grid = CPU()` or `arch_or_grid = GPU()`, data is interpolated
                  on-the-fly when the forcing tendency is computed. Default: CPU().

Keyword Arguments
=================

- `rate`: The restoring rate, i.e., the inverse of the restoring timescale (in s⁻¹).

- `mask`: The mask value. Can be a function of `(x, y, z, time)`, an array, or a number.

- `time_indices_in_memory`: The number of time indices to keep in memory. The number is chosen based on
                            a trade-off between increased performance (more indices in memory) and reduced
                            memory footprint (fewer indices in memory). Default: 2.

- `time_indexing`: The time indexing scheme for the field time series. Default: `Cyclical()`.

- `inpainting`: inpainting algorithm, see [`inpaint_mask!`](@ref). Default: `NearestNeighborInpainting(Inf)`.

- `cache_inpainted_data`: If `true`, the data is cached to disk after inpainting for later retrieving.
                          Default: `true`.
"""
function DatasetRestoring(metadata::Metadata,
                          arch_or_grid = CPU();
                          rate,
                          mask = 1,
                          time_indices_in_memory = default_time_indices_in_memory(metadata),
                          time_indexing = Cyclical(),
                          inpainting = NearestNeighborInpainting(Inf),
                          cache_inpainted_data = true)

    download_dataset(metadata)

    fts = FieldTimeSeries(metadata, arch_or_grid;
                          time_indices_in_memory,
                          time_indexing,
                          inpainting,
                          cache_inpainted_data)

    # Grab the correct Oceananigans field to restore
    variable_name = metadata.name
    field_name = oceananigans_fieldnames[variable_name]

    # If we pass the grid we do not need to interpolate
    # so we can save parameter space by setting the native grid to nothing
    on_native_grid = arch_or_grid isa AbstractArchitecture
    maybe_native_grid = on_native_grid ? fts.grid : nothing

    return DatasetRestoring(fts, maybe_native_grid, mask, field_name, rate)
end

function Base.show(io::IO, dsr::DatasetRestoring)
    print(io, "DatasetRestoring:", '\n',
              "├── variable_name: ", summary(dsr.variable_name), '\n',
              "├── rate: ", dsr.rate, '\n',
              "├── field_time_series: ", summary(dsr.field_time_series), '\n',
              "│    ├── dataset: ", summary(dsr.field_time_series.backend.metadata.dataset), '\n',
              "│    ├── dates: ", dsr.field_time_series.backend.metadata.dates, '\n',
              "│    ├── time_indexing: ", summary(dsr.field_time_series.time_indexing), '\n',
              "│    └── dir: ", dsr.field_time_series.backend.metadata.dir, '\n',
              "├── mask: ", summary(dsr.mask), '\n',
              "└── native_grid: ", summary(dsr.native_grid))
end

regularize_forcing(forcing::DatasetRestoring, field, field_name, model_field_names) = forcing

#####
##### Masks for restoring
#####

struct LinearlyTaperedPolarMask{N, S, Z}
    northern :: N
    southern :: S
    z :: Z
end

"""
    LinearlyTaperedPolarMask(; northern = (70,   75),
                               southern = (-75, -70),
                               z = (-20, 0))

Build a mask that is linearly tapered in latitude between the northern and southern edges.
The mask is constant in depth between the z and equals zero everywhere else.
The mask is limited to lie between (0, 1).
The mask has the following functional form:

```julia
n = 1 / (northern[2] - northern[1]) * (φ - northern[1])
s = 1 / (southern[1] - southern[2]) * (φ - southern[2])

valid_depth = (z[1] < z < z[2])

mask = valid_depth ? clamp(max(n, s), 0, 1) : 0
```
"""
function LinearlyTaperedPolarMask(; northern = (70,   75),
                                    southern = (-75, -70),
                                    z = (-20, 0))

    northern[1] > northern[2]  && throw(ArgumentError("Northern latitude range is invalid, northern[1] > northern[2]."))
    southern[1] > southern[2]  && throw(ArgumentError("Southern latitude range is invalid, southern[1] > southern[2]."))
    z[1] > z[2]                && throw(ArgumentError("Depth range is invalid, z[1] > z[2]."))

    return LinearlyTaperedPolarMask(northern, southern, z)
end

@inline function (mask::LinearlyTaperedPolarMask)(φ, z)
    # The mask is active only between `mask.z[1]` and `mask.z[2]`
    @inbounds begin
        northern_ramp = 1 / (mask.northern[2] - mask.northern[1]) * (φ - mask.northern[1])
        southern_ramp = 1 / (mask.southern[1] - mask.southern[2]) * (φ - mask.southern[2])
        within_depth_range = @inbounds (mask.z[1] < z < mask.z[2])
    end

    # Intersect the different mask components
    mask_value = max(northern_ramp, southern_ramp) * within_depth_range

    # Clamp the mask between 0 and 1
    mask_value = clamp(mask_value, 0, 1)

    return mask_value
end

@inline function stateindex(mask::LinearlyTaperedPolarMask, i, j, k, grid, time, loc)
    LX, LY, LZ = loc
    λ, φ, z = node(i, j, k, grid, LX(), LY(), LZ())
    return mask(φ, z)
end