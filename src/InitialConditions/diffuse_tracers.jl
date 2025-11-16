
scale_to_diffusivity(vertical_scale) = vertical_scale^2
scale_to_diffusivity(vertical_scale::Function) = (x, y, z, t) -> vertical_scale(x, y, z, t)^2

@kernel function _apply_tracer_mask!(tracer, initial_tracer, mask::AbstractArray)
    i, j, k = @index(Global, NTuple)
    @inbounds tracer[i, j, k] = ifelse(mask[i, j, k] == 1, initial_tracer[i, j, k], tracer[i, j, k])
end

@kernel function _apply_tracer_mask!(tracer, initial_tracer, mask::Function)
    i, j, k = @index(Global, NTuple)
    @inbounds tracer[i, j, k] = ifelse(mask(i, j, k, grid) == 1, initial_tracer[i, j, k], tracer[i, j, k])
end

@inline not_an_immersed_cell(i, j, k, grid) = !immersed_cell(i, j, k, grid)

function diffuse_tracers(initial_tracers;
                         horizontal_scale = 0,
                         vertical_scale = 0,
                         fractional_time_step = 2e-2,
                         mask = not_an_immersed_cell)

    # Remake the grid
    grid = initial_tracers[1].grid

    # Remake tracers without boundary conditions
    tracers = NamedTuple(name => CenterField(grid) for name in keys(initial_tracers))

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

    set!(smoothing_model, (tracer[name] = initial_tracer[name] for name in keys(initial_tracers))...)

    Nt = ceil(Int, 1 / Δt)
    @info string("Smoothing tracers ", keys(tracers), " with ", Nt, " time steps")
    smoothing_simulation = Simulation(smoothing_model; Δt, stop_time=1)

    # Remove NaN checker
    pop!(smoothing_simulation.callbacks, :nan_checker)

    # Restore values to default in the masked region
    function restore_values(sim)
        grid    = sim.model.grid
        tracers = sim.model.tracers
        for (tracer, initial_tracer) in zip(tracers, initial_tracers)
            launch!(architecture(grid), grid, :xyz, _apply_tracer_mask!, tracer, initial_tracer, mask)
        end
    end

    # Add restoring to initial values
    simulation.callbacks[:restoring] = Callback(restore_values, IterationInterval(1))

    run!(smoothing_simulation)

    return smoothing_model.tracers
end
