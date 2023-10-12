module InitialConditions

export continue_downwards!

using Oceananigans
using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior, CenterField
using Oceananigans.Architectures: architecture, device, GPU

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

function continue_downards!(field)
    arch = architecture(field)
    grid = field.grid
    loc = instantiated_location(field)
    launch!(arch, grid, :xy, _continue_downwards!, field, loc, grid)
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

scale_to_diffusivity(vertical_scale) = vertical_scale^2
scale_to_diffusivity(vertical_scale::Function) = (x, y, z, t) -> vertical_scale(x, y, z, t)^2

function diffuse_tracers!(grid;
                          tracers,
                          horizontal_scale = 0,
                          vertical_scale = 0,
                          fractional_time_step = 2e-2)

    # Remake tracers without boundary conditions
    tracers = NamedTuple(name => CenterField(grid; data=tracers[name].data) for name in keys(tracers))

    # Horizontal diffusivities that mix up to t ∼ ℓ² / κ ∼ 1
    κh = horizontal_scale^2
    κz = scale_to_diffusivity(vertical_scale)

    # Determine stable time-step
    Nx, Ny, Nz = size(grid)
    ϵ = fractional_time_step
    Az = minimum(grid.Azᶜᶜᵃ[1:Ny])

    if κh == 0
        Δt = fractional_time_step
    else
        Δt = ϵ * Az / κh
    end
    @show Nt = ceil(Int, 1 / Δt)

    vitd = VerticallyImplicitTimeDiscretization()
    vertical_smoothing = VerticalScalarDiffusivity(vitd, κ=κz)
    horizontal_smoothing = HorizontalScalarDiffusivity(κ=κh)

    smoothing_model = HydrostaticFreeSurfaceModel(; grid, tracers,
                                                  velocities = PrescribedVelocityFields(),
                                                  momentum_advection = nothing,
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

