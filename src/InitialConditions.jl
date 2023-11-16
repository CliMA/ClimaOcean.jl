module InitialConditions

export continue_downwards!

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields: OneField
using Oceananigans.Grids: peripheral_node
using Oceananigans.Utils: launch!
using Oceananigans.Fields: instantiated_location, interior, CenterField
using Oceananigans.Architectures: architecture, device, GPU

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

# Maybe we can remove this propagate field in lieu of a diffusion, 
# Still we'll need to do this a couple of steps on the original grid
@kernel function _propagate_field!(field, tmp_field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        nw = field[i - 1, j, k]
        ns = field[i, j - 1, k]
        ne = field[i + 1, j, k]
        nn = field[i, j + 1, k]
        nb = (nw, ne, nn, ns)

        counter = 0
        cumsum  = 0.0 

        @unroll for n in nb
            counter += ifelse(isnan(n), 0, 1)
            cumsum  += ifelse(isnan(n), 0, n)
        end

        tmp_field[i, j, k] = ifelse(cumsum == 0, NaN, cumsum / counter)
    end
end

@kernel function _substitute_values!(field, tmp_field)
    i, j, k = @index(Global, NTuple)
    @inbounds substitute = isnan(field[i, j, k])
    @inbounds field[i, j, k] = ifelse(substitute, tmp_field[i, j, k], field[i, j, k])
end

@kernel function _nans_outside_mask!(field, mask)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = ifelse(mask[i, j, k] == 0, NaN, field[i, j, k])
end

propagate_horizontally!(field, ::Nothing; kw...) = nothing

function propagate_horizontally!(field, mask; max_iter = Inf) 
    iter  = 0
    grid  = field.grid
    arch  = architecture(grid)
    
    launch!(arch, grid, :xyz, _nans_outside_mask!, field, mask)
    fill_halo_regions!(field)

    tmp_field = deepcopy(field)

    while isnan(sum(interior(field))) && iter < max_iter
        launch!(arch, grid, :xyz, _propagate_field!,   field, tmp_field)
        launch!(arch, grid, :xyz, _substitute_values!, field, tmp_field)
        iter += 1
        @debug "propagate pass $iter with sum $(sum(parent(field)))"
    end

    return nothing
end

continue_downwards!(field, ::Nothing) = nothing

function continue_downwards!(field, mask)
    arch = architecture(field)
    grid = field.grid
    launch!(arch, grid, :xy, _continue_downwards!, field, grid, mask)
    return nothing
end

@kernel function _continue_downwards!(field, grid, mask)
    i, j = @index(Global, NTuple)

    Nz = grid.Nz

    @unroll for k = Nz-1 : -1 : 1
        @inbounds fill_from_above = mask[i, j, k] == 0
        @inbounds field[i, j, k] = ifelse(fill_from_above, field[i, j, k+1], field[i, j, k])
    end
end

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

function adjust_tracers!(tracers; mask = nothing, max_iter = Inf)

    for (name, tracer) in zip(keys(tracers), tracers)
        @info "extending tracer $name"
        continue_downwards!(tracer, mask)
        propagate_horizontally!(tracer, mask; max_iter)
    end

    # Do we need this? Probably on the final grid
    # tracers = diffuse_tracers(tracers;
    #                           horizontal_scale,
    #                           vertical_scale,
    #                           fractional_time_step,
    #                           mask)

    return tracers
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

# TODO: move all the folowing to Oceananigans!

using Oceananigans.Fields: regrid!
using Oceananigans.Grids: cpu_face_constructor_x, 
                          cpu_face_constructor_y, 
                          cpu_face_constructor_z,
                          topology

# Should we move this to grids??
construct_grid(::Type{<:RectilinearGrid}, arch, size, extent, topology) = 
    RectilinearGrid(arch; size, x = extent[1], y = extent[2], z = extent[2], topology)

construct_grid(::Type{<:LatitudeLongitudeGrid}, arch, size, extent, topology) = 
    LatitudeLongitudeGrid(arch; size, longitude = extent[1], latitude = extent[2], z = extent[3], topology)

# Extend this to regrid automatically
function three_dimensional_regrid!(a, b)
    target_grid = a.grid isa ImmersedBoundaryGrid ? a.grid.underlying_grid : a.grid
    source_grid = b.grid isa ImmersedBoundaryGrid ? b.grid.underlying_grid : b.grid 

    topo = topology(target_grid)
    arch = architecture(target_grid)
    
    target_y = yt = cpu_face_constructor_y(target_grid)
    target_z = zt = cpu_face_constructor_z(target_grid)

    target_size = Nt = size(target_grid)

    source_x = xs = cpu_face_constructor_x(source_grid)
    source_y = ys = cpu_face_constructor_y(source_grid)

    source_size = Ns = size(source_grid)

    # Start by regridding in z
    @info "Regridding in z"
    zgrid   = construct_grid(typeof(target_grid), arch, (Ns[1], Ns[2], Nt[3]), (xs, ys, zt), topo)
    field_z = Field(location(b), zgrid)
    regrid!(field_z, zgrid, source_grid, b)

    # regrid in y 
    @info "Regridding in y"
    ygrid   = construct_grid(typeof(target_grid), arch, (Ns[1], Nt[2], Nt[3]), (xs, yt, zt), topo)
    field_y = Field(location(b), ygrid);
    regrid!(field_y, ygrid, zgrid, field_z);

    # Finally regrid in x
    @info "Regridding in x"
    regrid!(a, target_grid, ygrid, field_y)
end

end # module

