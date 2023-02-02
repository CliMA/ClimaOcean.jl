module DataWrangling

export continue_downwards!

using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior
using Oceananigans.Architectures: architecture, device_event, device, GPU

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

using ImageInpainting

function continue_downards!(field)
    arch = architecture(field)
    grid = field.grid
    loc = instantiated_location(field)

    event = launch!(arch, grid, :xy, _continue_downwards!, field, loc, grid; dependencies=device_event(arch))

    wait(device(arch), event)

    return nothing
end

@kernel function _continue_downwards!(field, (LX, LY, LZ), grid)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    active_surface = !peripheral_node(i, j, Nz, grid, LX, LY, LZ)

    @unroll for k = Nz-1 : -1 : 1
        fill_from_above = active_surface & peripheral_node(i, j, k, grid, LX, LY, LZ)
        @inbounds field[i, j, k] = ifelse(fill_from_above, field[i, j, k+1], field[i, j, k])
    end
end

function inpaint_horizontally!(field; algorithm=Criminisi(11, 11))
    arch = architecture(field)
    grid = field.grid
    loc = instantiated_location(field)
    Nx, Ny, Nz = size(grid)

    # Could transfer to CPU, then transfer back.
    arch isa GPU && error("Inpainting on the GPU is not supported yet!")

    for k = 1:Nz
        mask = [peripheral_node(i, j, k, grid, loc...) for i = 1:Nx, j = 1:Ny]
        interior(field, :, :, k) .= inpaint(interior(field, :, :, k), mask, algorithm) 
    end

    return nothing
end

end # module

