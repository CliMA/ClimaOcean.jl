using Oceananigans
using Oceananigans.Grids: halo_size
using Oceananigans.Fields: location
using JLD2

const Nx = 256
const Ny = 128
const Nz = 1
const nx = 256 ÷ 2
const ny = 128 ÷ 2

grid = TripolarGrid(size=(Nx, Ny, Nz))
Hx, Hy, Hz = halo_size(grid)

function read_bathymetry(prefix)
    bottom_height = zeros(Nx, Ny)

        for xrank in 0:1, yrank in 0:1
            rank = yrank + 2 * xrank
            file   = jldopen(prefix * "_$(rank).jld2")
            irange = nx * xrank + 1 : nx * (xrank + 1)
            jrange = ny * yrank + 1 : ny * (yrank + 1)

            data   = file["serialized/grid"].immersed_boundary.bottom_height[Hx+1:nx+Hx, Hy+1:ny+Hy, 1]
            bottom_height[irange, jrange] .= data
            close(file)
    end

    return bottom_height
end

bottom_height = read_bathymetry("surface_fields")

grid  = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

file0 = jldopen("surface_fields_0.jld2")
iters = keys(file0["timeseries/t"])
times = Float64[file0["timeseries/t/$(iter)"] for iter in iters]
close(file0)

utmp = FieldTimeSeries{Face,   Center, Nothing}(grid, times; backend=OnDisk(), path="surface_fields.jld2", name="u")
vtmp = FieldTimeSeries{Center, Face,   Nothing}(grid, times; backend=OnDisk(), path="surface_fields.jld2", name="v")
Ttmp = FieldTimeSeries{Center, Center, Nothing}(grid, times; backend=OnDisk(), path="surface_fields.jld2", name="T")
Stmp = FieldTimeSeries{Center, Center, Nothing}(grid, times; backend=OnDisk(), path="surface_fields.jld2", name="S")
etmp = FieldTimeSeries{Center, Center, Nothing}(grid, times; backend=OnDisk(), path="surface_fields.jld2", name="e")

function set_distributed_field_time_series!(fts, prefix)
    field = Field{location(fts)...}(grid)
    Ny = size(fts, 2)
    for (idx, iter) in enumerate(iters)
        @info "doing iter $idx of $(length(iters))"
        for xrank in 0:1, yrank in 0:1
            rank = yrank + 2 * xrank
            file   = jldopen(prefix * "_$(rank).jld2")
            irange = nx * xrank + 1 : nx * (xrank + 1)
            jrange = ny * yrank + 1 : ny * (yrank + 1)
            data   = file["timeseries/$(fts.name)/$(iter)"][Hx+1:nx+Hx, Hy+1:ny+Hy, 1]
            
            interior(field, irange, jrange, 1) .= data
            close(file)
        end

        set!(fts, field, idx)
    end
end

set_distributed_field_time_series!(utmp, "surface_fields")
set_distributed_field_time_series!(vtmp, "surface_fields")
set_distributed_field_time_series!(Ttmp, "surface_fields")
set_distributed_field_time_series!(Stmp, "surface_fields")
set_distributed_field_time_series!(etmp, "surface_fields")