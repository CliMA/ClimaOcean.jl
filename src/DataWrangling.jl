module DataWrangling

export continue_downwards!

using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior
using Oceananigans.Architectures: architecture, device_event, device, GPU

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

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

#=
using ImageInpainting

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
=#

function diffuse_tracers!(grid;
                          tracers,
                          horizontal_scale = 0,
                          vertical_scale = 0,
                          fractional_time_step = 2e-2)

    # Horizontal diffusivities that mix up to t ∼ ℓ² / κ ∼ 1
    κh = horizontal_scale^2
    κz = vertical_scale^2

    # Determine stable time-step
    grid = simulation.model.grid
    Nx, Ny, Nz = size(grid)
    ϵ = fractional_time_step
    Az = minimum(grid.Azᶜᶜᵃ[1:Ny])
    Δt = ϵ * Az / κh
    @show Nt = ceil(Int, 1 / Δt)

    vitd = VerticallyImplicitTimeDiscretization()
    vertical_smoothing = VerticalScalarDiffusivity(vitd, κ=κz)
    horizontal_smoothing = HorizontalScalarDiffusivity(κ=κh)

    smoothing_model = HydrostaticFreeSurfaceModel(; grid, tracers,
                                                  velocities = PrescribedVelocityFields(),
                                                  tracer_advection = nothing,
                                                  buoyancy = nothing,
                                                  closure = (horizontal_smoothing, vertical_smoothing))

    Nt = ceil(Int, 1 / Δt)
    @info string("Smoothing tracers ", keys(tracers), " with ", Nt, " time steps")
    smoothing_simulation = Simulation(smoothing_model; Δt, stop_time=1)
    pop!(smoothing_simulation.callbacks, :nan_checker) 
    run!(smoothing_simulation)

    return nothing
end

end # module

